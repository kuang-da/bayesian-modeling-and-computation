# Bayesian Modeling and Computation

## Course Topics

- [ ] 01. Introduction to Bayesian Probability Modeling (Ch.1,2,3)(Lec. 1-4) 
- [ ] 02. Regression Models from the Bayesian Perspective (Ch. 14,15)(Lec. 5,6) 
- [ ] 03. Hierarchical and Mixture Models (Ch. 5,22)(Lec. 7,8) 
- [ ] 04. Optimization for Model Estimation: EM and Variational Inference (Ch. 13)(Lec. 9) 
- [ ] 05. Markov Chain Monte Carlo for Model Estimation (Ch. 10,11)(Lec. 10, 11)
- [ ] 06. Recent Advances in Monte Carlo Simulation (Ch. 12)
- [ ] 07. Model Checking (Ch. 6,7)
- [ ] 08. General Linear Models (Ch. 16)
- [ ] 09. Hidden Markov Models
- [ ] 10. Dynamic Linear Models
- [ ] 11. Non-parametric Bayesian models (Ch. 23)
- [ ] 12. Gaussian Processes (Ch. 21)
- [ ] 13. Bayesian Tree Models

## Deployment

Tagging a commit will trigger the GitHub Action to generate a new release of PDF files.

To push a new commit with a tag,

```bash
git commit -m "prepare for v1.0.0 release"
git tag v1.0.0
git push origin main --tags
```
To tag an exist commit

```bash
git pull
git tag v0.0.3 44f171b3b29452269325c0363df0dc066c452154
git push origin main --tags
```
