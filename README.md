[![Crates.io](https://img.shields.io/crates/d/scarus.svg)](https://crates.io/crates/scarus)
[![Crates.io](https://img.shields.io/crates/v/scarus.svg)](https://crates.io/crates/scarus)
[![Crates.io](https://img.shields.io/crates/l/scarus.svg)](https://crates.io/crates/scarus)
[![CI](https://github.com/bihealth/scarus/actions/workflows/rust.yml/badge.svg)](https://github.com/bihealth/scarus/actions/workflows/rust.yml)
[![codecov](https://codecov.io/gh/bihealth/scarus/branch/main/graph/badge.svg?token=aZchhLWdzt)](https://codecov.io/gh/bihealth/scarus)

# Scarus

<img src="https://raw.githubusercontent.com/bihealth/scarus/main/utils/dna-fish.jpeg" width="256px" height="256px" align="right">

Scarus (Simple Code for ACMG Rule Utilisation and Scrutiny) is a software package for the autoamted evaluation of ACMG criteria for structural variants.

Also, here are some facts on Scarus/Parrotfish:

#1

> Scarus is a genus of parrotfishes.
>
> [...]
>
> Plutarch had written that scarus fish "swim together in shoals and ingeniously and heroically free each other when caught in a net."
> The scarus thus "denoted reciprocal assistance in the fight for survival."
>
> -- ["Scarus" -- Wikipedia](https://en.wikipedia.org/wiki/Scarus)

#2

> Their feeding activity is important for the production and distribution of coral sands in the reef biome, and can prevent algal overgrowth of the reef structure.
> The teeth grow continuously, replacing material worn away by feeding
>
> -- ["Parrotfish" -- Wikipedia](https://en.wikipedia.org/wiki/Parrotfish)

#3

> Corallivores are an important group of reef organism because they can influence coral abundance, distribution, and community structure.
>
> -- ["Corallivore" -- Wikipedia](https://en.wikipedia.org/wiki/Corallivore)

## Developer Documentation

The contents of this section are only relevant to developers of Scarus itself.

### Managing Project with Terraform

```
# export GITHUB_OWNER=bihealth
# export GITHUB_TOKEN=ghp_TOKEN

# cd utils/terraform
# terraform init
# terraform import github_repository.scarus scarus
# terraform validate
# terraform fmt
# terraform plan
# terraform apply
```
