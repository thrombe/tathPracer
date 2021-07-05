package main

import (
    "fmt"
    "time"
    "math"
    // "math/rand"
)

func shhh(vals ...interface{}) {
    for _, val := range vals {
        _ = val
    }
}

func main() {
    start := time.Now()
    defer func() {fmt.Println(time.Now().Sub(start))}()
    engine()
}

func camInit() *camera {
    cam := camera{}
    cam.height = 480
    cam.width = 720
    cam.fov = math.Pi/3
    cot := 1/math.Atan(cam.fov/2)
    cam.scrDist = (float64(cam.width)/2)*cot
    return &cam
}

type camera struct {
    height, width int
    fov, scrDist float64
    camPos, camDir [][]float64
}