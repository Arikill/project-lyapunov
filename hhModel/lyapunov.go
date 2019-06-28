package hhmodel

import (
	"fmt"
	"math"
)

func computeFutureStates(v float64) [4]float64 {
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

	dn := computeGateRate(n, nAlpha, nBeta)
	dm := computeGateRate(m, mAlpha, mBeta)
	dh := computeGateRate(h, hAlpha, hBeta)

	dv := gradientOfMembranePotential(v)
	futureStates := [4]float64{
		dv, dm, dh, dn,
	}

	return futureStates
}

func computeFutureStatesJacobian(states [4]float64) [4][4]float64 {
	var (
		dvBydv float64
		dvBydm float64
		dvBydh float64
		dvBydn float64

		dmBydv float64
		dmBydm float64
		dmBydh float64
		dmBydn float64

		dhBydv float64
		dhBydm float64
		dhBydh float64
		dhBydn float64

		dnBydv float64
		dnBydm float64
		dnBydh float64
		dnBydn float64
	)

	// First state (V) jacobian:
	dvBydv = -36*math.Pow(states[3], 4) - 120*states[2]*math.Pow(states[1], 3) - 3/10
	dvBydm = -144 * math.Pow(states[1], 3) * (states[0] + 77)
	dvBydh = -120 * math.Pow(states[3], 3) * (states[0] - 50)
	dvBydn = -360 * math.Pow(states[3], 2) * states[2] * (states[0] - 50)

	// Second state (m) jacobian:
	dmBydv = (states[3]*math.Exp(-states[0]/80-13/16))/640 + (states[1]-1)/(100*(math.Exp(-states[0]/10-11/2)-1)) + (math.Exp(-states[0]/10-11/2)*(states[0]/100+11/20)*(states[1]-1))/(10*math.Pow(math.Exp(-states[0]/10-11/2)-1, 2))
	dmBydm = (states[0]/100+11/20)/(math.Exp(-states[0]/10-11/2)-1) - math.Exp(-states[0]/80-13/16)/8
	dmBydh = 0
	dmBydn = 0

	// Third state (h) jacobian:
	dhBydv = (7*math.Exp(-states[0]/20-13/4)*(states[2]-1))/2000 - (states[2]*math.Exp(-states[0]/10-7/2))/(10*math.Pow(math.Exp(-states[0]/10-7/2)+1, 2))
	dhBydm = 0
	dhBydh = -(7*math.Exp(-states[0]/20-13/4))/100 - 1/(math.Exp(-states[0]/10-7/2)+1)
	dhBydn = 0

	// Fourth states (n) jacobian:
	dnBydv = (2*states[3]*math.Exp(-states[0]/18-65/18))/9 + (states[3]-1)/(10*(math.Exp(-states[0]/10-4)-1)) + (math.Exp(-states[0]/10-4)*(states[0]/10+4)*(states[3]-1))/(10*math.Pow(math.Exp(-states[0]/10-4)-1, 2))
	dnBydm = 0
	dnBydh = 0
	dnBydn = (states[0]/10+4)/(math.Exp(-states[0]/10-4)-1) - 4*math.Exp(-states[0]/18-65/18)

	result := [4][4]float64{
		{dvBydv, dvBydm, dvBydh, dvBydn},
		{dmBydv, dmBydm, dmBydh, dmBydn},
		{dhBydv, dhBydm, dhBydh, dhBydn},
		{dnBydv, dnBydm, dnBydh, dnBydn},
	}

	fmt.Println(result)

	return result
}
