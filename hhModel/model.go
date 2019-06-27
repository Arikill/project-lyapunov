package hhmodel

import "math"

const (
	// Er Resting Potential
	Er float64 = -54.387
	// ENa Sodium Reversal Potential
	ENa float64 = 50
	// EK Potassium Reversal Potential
	EK float64 = -77
	// Cm Membrane Capacitance
	Cm float64 = 1
	// gl leakage conductance
	gl float64 = 0.3
	// gNa peak sodium conductance
	gNa float64 = 120
	// gK peak potassium conductance
	gK float64 = 36
)

func computeAlpha(gateType string, v float64) float64 {
	if gateType == "n" {
		// alpha(n) is the number of times gate n opens per sec from a shut state
		return (-(0.01*v + 0.55) / (math.Exp(-v/10-11/2) - 1))
	} else if gateType == "m" {
		// alpha(m) is the number of times gate m opens per sec from a shut state
		return (0.1 * (v + 40) / (1 - math.Exp(-(v+40)/10)))
	} else if gateType == "h" {
		// alpha(h) is the number of times gate h opens per sec from a shut state
		return (0.07 * math.Exp(-(v+65)/20))
	}
	return 0
}

func computeBeta(gateType string, v float64) float64 {
	if gateType == "n" {
		// beta(n) is the number of times gate n shuts per sec from an open state
		return (0.125 * math.Exp(-(v+65)/80))
	} else if gateType == "m" {
		// beta(m) is the number of times gate m shuts per sec from an open state
		return (4 * math.Exp(-(v+65)/18))
	} else if gateType == "h" {
		// beta(h) is the number of times gate h shuts per sec from an open state
		return (1 / (1 + math.Exp(-(v+35)/10)))
	}
	return 0
}

// computeGateValue computes the probability that the gate is open
func computeGateValue(alpha float64, beta float64) float64 {
	return (alpha / (alpha + beta))
}

// computeGateRate computes the rate at which a gate changes its probability of opening.
func computeGateRate(gateValue float64, alpha float64, beta float64) float64 {
	return (alpha*(1-gateValue) + beta*gateValue)
}

// gradientOfMembranePotential computes dv/dt
func gradientOfMembranePotential(v float64) float64 {
	// Compute gate parameters
	nAlpha := computeAlpha("n", v)
	mAlpha := computeAlpha("m", v)
	hAlpha := computeAlpha("h", v)
	nBeta := computeBeta("n", v)
	mBeta := computeBeta("m", v)
	hBeta := computeBeta("h", v)

	// Compute gate values
	n := computeGateValue(nAlpha, nBeta)
	m := computeGateValue(mAlpha, mBeta)
	h := computeGateValue(hAlpha, hBeta)

	return ((1 / Cm) * ((-gNa * math.Pow(m, 3) * h * (v - ENa)) - (gK * math.Pow(n, 4) * (v - EK)) - (gl * (v - Er))))
}
