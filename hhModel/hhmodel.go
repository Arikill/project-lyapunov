package hhmodel

// ComputeEquilibrium a wrapper computes equilibrium points for hhmodel.
func ComputeEquilibrium() []float64 {
	return equilibriumPointsGradientDecent()
}

// ComputeLyapunovStability a wrapper that computes lyapunov stability of hhmodel.
func ComputeLyapunovStability() {

}
