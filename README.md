
# MRSS - Multi Resolution Subsampled Saliency

Code for the paper "Saliency detection for large-scale mesh decimation", C&G 2023
https://www.sciencedirect.com/science/article/pii/S0097849323000134?via%3Dihub

Implements the method described in the paper and a decimation algorithm that uses the saliency information. Saliency maps can also be used on ZBrush as the "polypaint" if you are working with a quadmesh. 

Feel free to use it following the provided license (authors would be happy if you notify us when using it).

![Alt text](/relative/path/to/img.jpg?raw=true "Optional Title")



## Building

We use Cmake (https://cmake.org/) to generate project files. Usage would be from a terminal at the root folder of the project:

```bash
  mkdir build
  cmake ..
```
    
For now it creates a project called "green" (old temporary name). Will be updated to MRSS soon. 
## References

We have used:
assimp 5.0.1, fmt 7.1.3, glew 2.1.0, glfw 3.3, glm 0.9.9.5 ,imgui 1.74

OpenMesh (https://gitlab.vci.rwth-aachen.de:9000/OpenMesh/OpenMesh)

We also thank Max Limper for providing code for his EG paper (https://diglib.eg.org/handle/10.2312/egsh20161003) which greatly helped with our development. 

We would like to call special thanks to Ayumi Kimura for her continuous support. We are also grateful to Square Enix and OLM Digital for their thought-provoking discussions on our work. We also thank Warren Butcher and other peers from CMIC, for their
support during the development of this publication. 
