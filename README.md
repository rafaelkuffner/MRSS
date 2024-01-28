
# MRSS - Multi Resolution Subsampled Saliency
![Results](/images/header.png)
Code for the paper ["Saliency detection for large-scale mesh decimation", C&G 2023](https://www.sciencedirect.com/science/article/pii/S0097849323000134?via%3Dihub)

Implements the method described in the paper and a decimation algorithm that uses the saliency information. Saliency maps can also be used on ZBrush as the "polypaint" if you are working with a quadmesh. 
![GUI](/images/GUI.png)


Feel free to use it following the provided license (authors would be happy if you notify us when using it).




## Building

We use [CMake](https://cmake.org/) to generate project files. Usage would be from a terminal at the root folder of the project:

```bash
  mkdir build
  cmake -B build -S .
```
    
For now it creates a project called "green" (old temporary name). ~Will~ May be updated to MRSS ~soon~ eventually.

## References

Libraries used:
- [OpenMesh](https://www.graphics.rwth-aachen.de/software/openmesh/) (with various modifications)
- [AssImp](https://github.com/assimp/assimp)
- [GLFW](https://www.glfw.org/)
- [GLEW](https://glew.sourceforge.net/)
- [GLM](https://github.com/g-truc/glm)
- [ImGui](https://github.com/ocornut/imgui/)
- [fmt](https://github.com/fmtlib/fmt/)
- [clipp](https://github.com/muellan/clipp)

We also thank Max Limper for providing code for his [EG paper](https://diglib.eg.org/handle/10.2312/egsh20161003) which greatly helped with our development. 

We would like to call special thanks to Ayumi Kimura for her continuous support. We are also grateful to Square Enix and OLM Digital for their thought-provoking discussions on our work. We also thank Warren Butcher and other peers from CMIC, for their
support during the development of this publication. 
