package main

import (
	"fmt"
	"math"
)

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

func main() {
	eqPt := equilibriumPointsGradientDecent()
	fmt.Println("The eq points are: ", eqPt)
	fmt.Println("The value at equilibrium is: ", gradientOfMembranePotential(eqPt[0]))
}

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

func computeGateValue(alpha float64, beta float64) float64 {
	// Gate value (probablity that the gate is open)
	return (alpha / (alpha + beta))
}

func computeGateRate(gateValue float64, alpha float64, beta float64) float64 {
	// Rate at which a gate changes its probability of opening.
	return (alpha*(1-gateValue) + beta*gateValue)
}

// gradientOfMembranePotential computes
func gradientOfMembranePotential(v float64) float64 {
	nAlpha := computeAlpha("n", v)
	mAlpha := computeAlpha("m", v)
	hAlpha := computeAlpha("h", v)

	nBeta := computeBeta("n", v)
	mBeta := computeBeta("m", v)
	hBeta := computeBeta("h", v)

	n := computeGateValue(nAlpha, nBeta)
	m := computeGateValue(mAlpha, mBeta)
	h := computeGateValue(hAlpha, hBeta)

	return ((1 / Cm) * ((-gNa * math.Pow(m, 3) * h * (v - ENa)) - (gK * math.Pow(n, 4) * (v - EK)) - (gl * (v - Er))))
}

func gradientOfGradientOfMembranePotential(v float64) float64 {
	dv := gradientOfMembranePotential(v)
	nAlpha := computeAlpha("n", v)
	mAlpha := computeAlpha("m", v)
	hAlpha := computeAlpha("h", v)
	nBeta := computeBeta("n", v)
	mBeta := computeBeta("m", v)
	hBeta := computeBeta("h", v)
	n := computeGateValue(nAlpha, nBeta)
	m := computeGateValue(mAlpha, mBeta)
	h := computeGateValue(hAlpha, hBeta)
	dn := computeGateRate(n, nAlpha, nBeta)
	dm := computeGateRate(m, mAlpha, mBeta)
	dh := computeGateRate(h, hAlpha, hBeta)
	gNaCoefficients := -(3 * math.Pow(m, 2) * dm * h * (v - ENa)) - (math.Pow(m, 3) * dh * (v - ENa)) - (math.Pow(m, 3) * h * dv)
	gKCoefficients := -(4 * math.Pow(n, 3) * dn * (v - EK)) - (math.Pow(n, 4) * dv)
	glCoefficients := -dv
	return ((1 / Cm) * (gNaCoefficients*gNa + gKCoefficients*gK + glCoefficients*gl))
}

func objectiveFunction(v float64) float64 {
	return (0.5 * math.Pow(gradientOfMembranePotential(v), 2))
}

func gradientOfObjectiveFunction(v float64) float64 {
	return (gradientOfGradientOfMembranePotential(v) * gradientOfMembranePotential(v))
}

func equilibriumPointsGradientDecent() []float64 {
	var (
		vCurrent          = Er
		costCurrent       float64
		vNext             float64
		costNext          float64
		learningRate      = 0.01
		precision         = 1E-10
		iteration         int
		equilibriumPoints = make([]float64, 4)
	)
	for true {
		costCurrent = gradientOfMembranePotential(vCurrent)
		vNext = vCurrent - learningRate*gradientOfMembranePotential(vCurrent)
		costNext = gradientOfMembranePotential(vNext)
		fmt.Println(iteration, ": From ", vCurrent, " to ", vNext, " @ ", learningRate, " of cost: ", costCurrent)
		learningRate = (vNext - vCurrent) / math.Abs(gradientOfObjectiveFunction(vNext)-gradientOfObjectiveFunction(vCurrent))
		if math.Abs(costNext-costCurrent) <= precision || math.IsNaN(vNext) || math.IsNaN(learningRate) {
			break
		}
		vCurrent = vNext
		iteration++
	}
	nAlpha := computeAlpha("n", vNext)
	mAlpha := computeAlpha("m", vNext)
	hAlpha := computeAlpha("h", vNext)
	nBeta := computeBeta("n", vNext)
	mBeta := computeBeta("m", vNext)
	hBeta := computeBeta("h", vNext)
	n := computeGateValue(nAlpha, nBeta)
	m := computeGateValue(mAlpha, mBeta)
	h := computeGateValue(hAlpha, hBeta)
	equilibriumPoints[0] = vNext
	equilibriumPoints[1] = m
	equilibriumPoints[2] = h
	equilibriumPoints[3] = n
	return equilibriumPoints
}
