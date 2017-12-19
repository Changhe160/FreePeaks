# Free Peaks
This is an open source project with the aim of providing a unified and fully tunable framework for constructing continuous optimization problems. This framework is easy to understand, highly configurable, feature enriched, and computationally efficient. It is easy for users to analyze weaknesses and strengths of an algorithm. And you can construct a problem with desired features and with your own component functions. 


The framework utilizes the k-d tree to partition the search space and sets a certain number of simple functions in each subspace. And it can be used for constructing global/multimodal optimization, dynamic single objective optimization, multi-objective optimization, and dynamic multi-objective optimization, which have the following features.

1. The framework is easy to understand and analyze in the following aspects: the position of an optimum, the shape of an optimum, the basin of attraction of an optimum, the local structure of an optimum, the shape of the Pareto optimal front (POF), the distribution of the Pareto optimal set (POS), and the dominating relationship between any solutions.
2. The framework is independently configurable in the following aspects: the dimensionality, the modality, the number of objectives, the position of an optimum, the shape of an optimum, the size of the basin of attraction of an optimum, the local/global structure, the shape of the POF, the distribution of the POS, the size of a countable POS, the inter-relationship between decision variables, the size of feasible areas, and domino convergence.
3. The framework is computationally efficient for fitness evaluations in comparison with existing problems.
4. The framework is able to show common characteristics of existing problems, e.g., the common shapes of the POF.
5. The framework is flexible and extendable, i.e., users are able to easily add new components to the framework.
6. The framework is able to easily manipulate any optimum with different transformations, e.g., rotation, shift, irregularities, and break of the symmetry of symmetric peaks.
7. The framework is compatible with existing problems. Existing problems can be integrated into the framework without changing their structures.

The detail of the framework is described in the paper below:
Changhe Li, Trung Thanh Nguyen, Sanyou Zeng, Ming Yang, Min Wu, "An Open Framework for Constructing Continuous Optimization Problems", to be appeared on IEEE Transactions on Cybernetics. 
  
Two example algorithms are available to show how to use the framework. One algorithm is the CMAES for global optimization and the other is the MOEA/D for multi-objective optimization. 
 
Free Peaks was implemented in C++11 language based on some libraries of Boost and several modern C++ design techniques. It can be run in both Windows and Linux environments with a single thread or multi-thread mode. The current version of Free Peaks is v0.1 developed with MS Visual Studio Community 2015.

This project is totally for research purpose. I would like to encourage you to take part in this project and let's together make the research life easy, comfortable, and efficient for you and for everyone. Please feel free to contact me by changhe.lw@gmail.com if you have any questions.