# about
this is a multithreaded cpu path tracer.
i gave up on the go implementation because of the speed of my implementation but the rust implimentation was bareable

# some images rendered using tathpracer
![spheres sunny](./images/image1629809159.jpeg)
![spheres night](./images/image1630064890.jpeg)
![glassy cube](./images/image1630845677.jpeg)
![voxel teapot](./images/image1631196071.jpeg)
![triangle glass teapot](./images/image1631545550.jpeg)
![voxel glass teapot](./images/image1631546482.jpeg)

# how to run
- clone repo
- cd to rust dir
- run ```cargo run --release```
- try different scenes by uncommenting/commenting lines src/scene.rs

# features
- loosely follows the theory from [ray tracing in one weekend](https://github.com/RayTracing/raytracing.github.io) with following extra features
- SVO implimentation + ray intersections
- different materials: metallic, diffuse, dielectric, glowing matt surfaces
- objects: sphere, plane, triangles, triangle octree, voxel octree
- import .obj triangle meshes
  - only triangle meshes supported 
  - colors not supported
- triangles are stored in SVO after importing
- simple conversion from triangle meshes to SVO