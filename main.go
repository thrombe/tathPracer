package main

import (
    "fmt"
    "time"
    "math"
)

func shhh(vals ...interface{}) {
    for _, val := range vals {
        _ = val
    }
}

func main() {
    start := time.Now()
    defer func() {fmt.Println(time.Now().Sub(start))}()
    tath()
}

func tath() {
    height, width := 480, 720
    img, set := newImg(width, height)
    fov := math.Pi/3
    cot := 1/math.Atan(fov/2)
    scrDist := (float64(width)/2)*cot
    
    // fwd := vector(0, 0, -1)
    rht := vector(1, 0, 0)
    up := vector(0, 1, 0)
    topLeft := vector(-float64(width)/2, float64(height)/2, -scrDist)
        
    camPos := vector(0, 0, 0)
    hitSP1 := sphereHitCreate(vector(0, 0, -5), 1)

    for y := 0; y < height; y++ {
        for x := 0; x < width; x++ {
            ray := nMatAdd(camPos, topLeft, matScalar(rht, float64(x)), matScalar(up, float64(-y)))
            bright, did := hitSP1(camPos, ray)
            if !did {
                backgroundPix(x, y, height, set)
                continue
            }
            // if y == 480/2 {fmt.Println(bright, intersectionPoint, t, x, y)}
            set(x, y, round(absVal(255*bright)), 0, 0)
        }
    }
    dumpImg(img)
}

func backgroundPix(x, y, height int, set func(int, int, int, int, int)) {
    t := 1 - float64(y)/float64(height)
    rad := 229*(1-t) + t*148
    grn := 240*(1-t) + t*191
    blu := 255*(1-t) + t*255
    set(x, y, int(rad), int(grn), int(blu))
}

func sphereHitCreate(center [][]float64, r float64) func([][]float64, [][]float64) (float64, bool) {
    return func(camPos, ray [][]float64) (float64, bool) {
        oc := matSub(camPos, center)
        negB := -vecDot(ray, oc)
        bSq := vecDot(ray, ray)
        Dby4 := negB*negB - bSq*(vecDot(oc, oc)-r*r)
        if Dby4 < 0 {return 0, false} // didnt hit
        t := (negB - math.Sqrt(Dby4))/bSq // no +ve sqrt(D) cuz we want min anyway
        if t < 0 {return 0, false} // ray hitting behind the camera
        intersectionPoint := matScalar(ray, t)
        intersectionNormal := matSub(intersectionPoint, center)
        intersectionNormal = vecUnit(intersectionNormal)
        // set(x, y, round(absVal(255*intersectionNormal[0][0])), round(absVal(255*intersectionNormal[1][0])), round(absVal(255*intersectionNormal[2][0])))
        bright := vecDot(intersectionNormal, vecUnit(vector(-1, -1, -2)))
        return bright, true
    }
}