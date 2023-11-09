[![Crates.io](https://img.shields.io/crates/d/scarus.svg)](https://crates.io/crates/scarus)
[![Crates.io](https://img.shields.io/crates/v/scarus.svg)](https://crates.io/crates/scarus)
[![Crates.io](https://img.shields.io/crates/l/scarus.svg)](https://crates.io/crates/scarus)
[![CI](https://github.com/bihealth/scarus/actions/workflows/rust.yml/badge.svg)](https://github.com/bihealth/scarus/actions/workflows/rust.yml)
[![codecov](https://codecov.io/gh/bihealth/scarus/branch/main/graph/badge.svg?token=aZchhLWdzt)](https://codecov.io/gh/bihealth/scarus)

# SCARUS - Simple Code for ACMG Rule Utilisation and Scrutiny

This repository contains code for the automated evaluation of ACMG criteria for sequence and structural variants.

## Managing Project with Terraform

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
