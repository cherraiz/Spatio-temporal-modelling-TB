library(INLA)

# counties_d <- shapefile with the polygons of the counties
# data_std <- main database with response and predictors
# data_wb10_std <- main database with response and predictors. Records with less than 10 wild boar tested excluded.

# Creating a counties' id for the matrixes. It must match in the main database
# Must be numeric for matching with the neighborhood matrix index
counties_d$county_id <- as.numeric(rownames(counties_d))

# Assign the same index in the database
data_std <- data_std %>% group_by(county) %>% merge(data.frame(counties_d)[c("county","county_id")]) 

# Centroid data for each county
counties_sp <- as(st_transform(counties_d, crs = 25830), "Spatial")
counties_sp$x <- coordinates(counties_sp)[,1]; counties_sp$y <- coordinates(counties_sp)[,2]
data_std <- data_std %>% group_by(county, county_id) %>% merge(counties_sp@data, all =T)


## 4.1.1. iCAR: Neighborhood matrix 
# Neighborhood list
nb <- poly2nb(counties_d, row.names = counties_d$county_id)

# Neighborhood matrix  
inla.nb <- as(nb2mat(nb, style = "B", zero.policy = TRUE), "Matrix") 

# Plot
plot(st_geometry(counties_d), col = "grey", main = "Neighborhood structure")
plot(nb, coordinates(as(counties_d, "Spatial")), col='red', lwd=1.2, add=TRUE, pch = 16)


## 4.1.2. SPDE: Simple Delaunay triangulation
mesh_sp <- inla.mesh.create(loc = cbind(data_std$x,data_std$y))
spde_sp <- inla.spde2.matern(mesh_sp)

# Plot
plot(mesh_sp)

# Mesh field in the database
data_std$mesh_sp <- mesh_sp$id$loc


## 4.1.3. CRDT: Constrained Refined Delaunay Triangulation
# For obtaining triangles as similar as possible in shape and size
# Maximum side size (max.edge) was calculated as 0.05 quantile of the side size of the counties considering them as triangles
# The minimum angle (min.angle) in 21, max. value guaranteeing the convergence of the algorithm (Krainski et al. 2018)
counties_re <- st_transform(counties_d, crs = 25830)
counties_re$area <- st_area(counties_re)
counties_re$triangle_side <- (sqrt(counties_re$area))/(sqrt(3)/4)
max.edge = quantile(counties_re$triangle_side, 0.05)
mesh_re <- inla.mesh.2d(loc = cbind(data_std$x, data_std$y), max.edge = max.edge, offset = c(150000,200000),min.angle = 21)
spde_re <- inla.spde2.matern(mesh_re)

# Plot
plot(mesh_re); plot(raster::aggregate(as(counties_re, "Spatial"), dissolve = TRUE), add=TRUE)

# Mesh field in the database
data_std$mesh_re <- mesh_re$id$loc


## 4.1.4. CRDT by Year
# For creating a spatio-temporal effect combining the structured spatial effect and the Year we need to create an inla.stack object with all data and the mesh
mesh_t <- inla.mesh.1d(loc = data_std$year)
spde_re <- inla.spde2.matern(mesh_re, alpha = 2)
A_re <- inla.spde.make.A(mesh = mesh_re, loc = cbind(data_wb10_std$x,data_wb10_std$y), group = data_wb10_std$year,
                         group.mesh = mesh_t)
index_re <- inla.spde.make.index(name = "spatial.effect", n.spde = spde_re$n.spde, n.group = mesh_t$m)
stack_re<- inla.stack(data=list(positive_herds = data_wb10_std$positive_herds, new_positive_herds = data_wb10_std$new_positive_herds,tested_herds = data_wb10_std$tested_herds), A=list(A_re, 1),   effects=list(c(index_re,list(Intercept=1)),list(data_wb10_std[,c("beef_herds","dairy_herds","bullfight_herds","HR_movements","LR_movements","HER_movements","ELR_movements","wild_boar_ab","red_deer_ab","TB_wild_boar")])), tag = "est")

                    