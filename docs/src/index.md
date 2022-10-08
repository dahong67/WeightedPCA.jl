```@meta
CurrentModule = WeightedPCA
```

# WeightedPCA: PCA for heterogeneous quality samples

Documentation for [WeightedPCA](https://github.com/dahong67/WeightedPCA.jl).

> 👋 *This package provides research code and work is ongoing.
> If you are interested in using it in your own research,
> **I'd love to hear from you and collaborate!**
> Feel free to write: [dahong67@wharton.upenn.edu](mailto:dahong67@wharton.upenn.edu)*

Please cite the following paper for this technique:

> David Hong, Fan Yang, Jeffrey A. Fessler, Laura Balzano.
> "Optimally Weighted PCA for High-Dimensional Heteroscedastic Data", 2022.
> [https://arxiv.org/abs/1810.12862](https://arxiv.org/abs/1810.12862)

In BibTeX form:
```bibtex
@Misc{hyfb2022owp,
   title  = "Optimally Weighted PCA for High-Dimensional Heteroscedastic Data", 
   author = "David Hong and Fan Yang and Jeffrey A. Fessler and Laura Balzano",
   year   = 2022,
   url    = "https://arxiv.org/abs/1810.12862",
}
```

## Docstrings

```@index
```

```@docs
WeightedPCA
wpca
ComputedWeights
UniformWeights
InverseVarianceWeights
OptimalWeights
```
