# Flow Field
 Visualisation of 3D flow fields using Julia. Intention of the project is for the user to enter a system of 3 ODEs for the flows to be visualised as an interactive 3D animation.
 
 Inspired heavily by [this project](anvaka.github.io/fieldplay) which visualises the phase portrait of 2D systems. I always found this project very useful during my studies, and wanted to extend it to work with an additional variable. Extending it in this way allows for the visualisation of interesting mathematical phenomena such as attractors.

 Initial work on this project was conducted using Python, however the available packages lacked the interactivity I wanted in the final application. In addition to this, the higher performance of Julia was more suitable for looped mathematical operations on large matrices.

 Within this repository includes a custom interface using the metaprogramming features of Julia. This allows an ODE to be fully specified using a simple `@ODE` macro.
 
 Future Work includes:
  - Implementation of ModelingToolkit.jl for more potential to integrate symbolic computation.
  - Additional options for interactivity, such as changing parameter values mid-visualisation.