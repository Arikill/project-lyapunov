package hhmodel

import "math"

// ComputeEquilibrium a wrapper computes equilibrium points for hhmodel.
func ComputeEquilibrium() []float64 {
	return equilibriumPointsGradientDecent()
}

// ComputeLyapunovStability a wrapper that computes lyapunov stability of hhmodel.
func ComputeLyapunovStability(v float64) {
	futureStates := computeFutureStates(v)
	jacobian := computeFutureStatesJacobian(futureStates)
	stability(int(math.Sqrt(float64(len(jacobian)))), int(math.Sqrt(float64(len(jacobian)))), jacobian)
	lipschitz(0, ENa)
}
