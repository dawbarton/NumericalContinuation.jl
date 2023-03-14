# NumericalContinuation

[![Build Status](https://github.com/dawbarton/NumericalContinuation.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/dawbarton/NumericalContinuation.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/dawbarton/NumericalContinuation.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/dawbarton/NumericalContinuation.jl)
[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet)](https://github.com/SciML/ColPrac)
[![SciML Code Style](https://img.shields.io/static/v1?label=code%20style&message=SciML&color=9558b2&labelColor=389826)](https://github.com/SciML/SciMLStyle)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

**If you are looking for a numerical continuation code in Julia, you are best off using [BifurcationKit.jl](https://github.com/bifurcationkit/BifurcationKit.jl) at the moment.**

This is a work-in-progress implementation of a numerical continuation code in the style of COCO (see Recipes for Continuation by Harry Dankowicz and Frank Schilder). The main difference to other numerical continuation approaches is the use of an extended state (and corresponding monitor functions) to provide considerable flexibility and ease of problem creation. Also, the notion of "staged construction" of problems is used; continuation problems can be built up in a hierarchical manner, e.g., a limit cycle problem can be composed with an equilibrium problem and a connecting orbit problem to track heteroclinic orbits (the individual parts are not all implemented yet!).

This is not a direct re-implementation of COCO in Julia; I am leveraging the Julia ecosystem and Julia idioms to avoid "re-implementing the wheel". Specifically, this is intended to work seemlessly with the [SciML/DiffEq ecosystem](https://sciml.ai/) and follows the same conventions wherever possible. Ultimately, my aim is to have numerical continuation tools that can be used in conjunction with hybrid machine learnt/physics-based models, and design decisions are made accordingly.

At the moment, I have implemented

- the staged-construction mechanism for continuation problems;
- optimised function generation (using `@generated` functions);
- monitor functions; and
- a simple Fourier-based collocation scheme for limit cycles.

Contributions are welcome even at this early stage. If you are interested in this project, please get in touch with [David Barton](mailto:david.barton@bristol.ac.uk).
