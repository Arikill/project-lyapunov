package main

import (
	"fmt"

	hhmodel "./hhmodel"
)

func main() {
	eqpts := hhmodel.ComputeEquilibrium()
	fmt.Println("V: ", eqpts[0])
	fmt.Println("m: ", eqpts[1])
	fmt.Println("h: ", eqpts[2])
	fmt.Println("n: ", eqpts[3])
	hhmodel.ComputeLyapunovStability(eqpts[0])
}
