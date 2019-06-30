package hhmodel

import (
	"fmt"
	"math"
)

// gradientOfGradientOfMembranePotential computes d2v/dt2
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
		// fmt.Println("Current cost: ", costCurrent)
		learningRate = (vNext - vCurrent) / math.Abs(gradientOfObjectiveFunction(vNext)-gradientOfObjectiveFunction(vCurrent))
		if math.Abs(costNext-costCurrent) <= precision || math.IsNaN(vNext) || math.IsNaN(learningRate) {
			fmt.Println("Converged! Minimum cost: ", costNext, " @ ", vNext)
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
