# MonteCarloPathTracer
A global illumination renderer using unbiased Monte Carlo path tracing.  
Based on <b>Assimp</b> and <b>glm</b>.  
Reference to <b>smallpt</b> and <b>Physically Based Rendering: From Theory to Implementation</b> and some other sites.  

## Input
- A model file and corresponding material file(you need to make sure they are no problem).

## Output
- A image.
- ![scene1](https://github.com/AmazingZhen/MonteCarloPathTracer/blob/master/MonteCarloPathTracer/resScene1/image_iterations_25000.png?raw=true)
- ![scene2](https://github.com/AmazingZhen/MonteCarloPathTracer/blob/master/MonteCarloPathTracer/resScene2/image_iterations_25000.png?raw=true)

## Algorithmic process
- Recursively ray tracing
  + Kd-tree for accelerating ray intersection
  + Monte-Carlo method for deciding next ray to trace
  + Different optical equation for final intensity calculation

## Areas for improvement
- Kd tree acceleration maybe exist improvement
- GPU based parallel ray tracing(very difficult)
- A better way to using Monte-Carlo method to sample, such as importance sample
