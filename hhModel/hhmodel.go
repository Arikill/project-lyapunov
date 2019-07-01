package hhmodel

import (
	"math"
)

// ComputeEquilibrium a wrapper computes equilibrium points for hhmodel.
func ComputeEquilibrium() []float64 {
	return equilibriumPointsGradientDecent()
}

// ComputeSystemStability a wrapper that computes lyapunov stability of hhmodel.
func ComputeSystemStability(v float64) bool {
	futureStates := computeFutureStates(v)
	jacobian := computeFutureStatesJacobian(futureStates)
	if stability(int(math.Sqrt(float64(len(jacobian)))), int(math.Sqrt(float64(len(jacobian)))), jacobian) && lipschitz(0, ENa) {
		return true
	}
	return false
}
