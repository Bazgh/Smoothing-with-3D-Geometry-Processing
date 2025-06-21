This project implements various mesh smoothing techniques to denoise and enhance 3D models while preserving important geometric features. The algorithms are tested on several mesh scans such as "Bunny" and "Max".

🧮 Implemented Methods
1. Explicit Smoothing
Uses the Laplacian operator to reduce high-frequency noise by updating vertex positions based on neighbors.

Uniform Smoothing:
Weights all neighboring vertices equally. Fast but may fail to preserve sharp features.

Cotangent Smoothing:
Weights neighbors based on cotangent values, allowing sharper edge preservation compared to uniform smoothing.

2. Implicit Smoothing
Solves a linear system involving the Laplacian for stable, large-step smoothing. This method:

Preserves surface convexity (analogous to the Gage-Hamilton-Grayson curvature flow theorem).

Effectively distinguishes noise from important geometric features.

Is computationally heavy and less effective on dense meshes.

3. Feature Enhancement
To prevent loss of important details in explicit methods:

A feature enhancement loop tracks large deviations between original and smoothed meshes.

Repeated smoothing + enhancement improves sharp feature preservation.

Optimal results achieved after ~3 iterations; more can introduce artifacts.

🔧 Parameters
time_step: Controls the smoothing aggressiveness.

iterations: Number of smoothing steps.

enhancement_coefficient: Controls the emphasis on sharp features.

Note: Over-tuning may lead to unrealistic meshes or loss of detail.

📉 Image Quality Impact
Higher-quality meshes (e.g., "nice max") yield better enhancement results. Poor-quality meshes may introduce artifacts due to noise misinterpretation.

🖼️ Visual Results
Each result image includes:

Uniformly smoothed & enhanced

Cotangent smoothed & enhanced

Original mesh
![Bunny Enhanced](uniform-cotan-original(bunny).png)
![Max Enhanced](uniform-cotan-original(nice max).png)
