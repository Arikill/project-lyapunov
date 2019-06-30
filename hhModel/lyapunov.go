package hhmodel

import (
	"fmt"
	"log"
	"math"

	"gonum.org/v1/gonum/mat"
)

func computeFutureStates(v float64) []float64 {
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

	futureStates := make([]float64, 4)
	futureStates[0] = dv
	futureStates[1] = dm
	futureStates[2] = dh
	futureStates[3] = dn

	return futureStates
}

func computeFutureStatesJacobian(states []float64) []float64 {
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

	// result := [16]float64{
	// 	dvBydv, dvBydm, dvBydh, dvBydn,
	// 	dmBydv, dmBydm, dmBydh, dmBydn,
	// 	dhBydv, dhBydm, dhBydh, dhBydn,
	// 	dnBydv, dnBydm, dnBydh, dnBydn,
	// }

	result := make([]float64, 16)
	result[0] = dvBydv
	result[1] = dvBydm
	result[2] = dvBydh
	result[3] = dvBydn

	result[4] = dmBydv
	result[5] = dmBydm
	result[6] = dmBydh
	result[7] = dmBydn

	result[8] = dhBydv
	result[9] = dhBydm
	result[10] = dhBydh
	result[11] = dhBydn

	result[12] = dnBydv
	result[13] = dnBydm
	result[14] = dnBydh
	result[15] = dnBydn

	fmt.Println(result)

	return result
}

func stability(rows int, cols int, slice []float64) []complex128 {
	// Generate a 6×6 matrix of random values.
	a := mat.NewDense(rows, cols, slice)
	var eig mat.Eigen
	ok := eig.Factorize(a, mat.EigenLeft)
	if !ok {
		log.Fatal("Eigendecomposition failed")
	}
	eigenValues := eig.Values(nil)
	fmt.Printf("Eigenvalues of A:\n%v\n", eigenValues)
	for _, value := range eigenValues {
		if real(value) > 0 {
			panic("system is unstable!")
		}
	}
	return eigenValues
}

func lipschitz(v1 float64, v2 float64) bool {
	dv1 := gradientOfMembranePotential(v1)
	dv2 := gradientOfMembranePotential(v2)
	fmt.Println(v1)
	fmt.Println(dv1)
	if math.Sqrt(math.Pow(dv1-dv2, 2)) <= math.Sqrt(math.Pow(v1-v2, 2)) {
		fmt.Println("Continuously differentiable!")
		return true
	}
	fmt.Println("Not Continuously differentiable!")
	return false
}
