package main

import (
    "fmt"
    "time"
    "math"
    "math/rand"
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

type camera struct {
    height, width int
    fov, scrDist float64
    camPos, camDir [][]float64
}

type pixel struct {
    x, y int
    r, g, b int
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

func tath() {
    cam := camInit()
    
    sem := make(chan struct{}, 8)
    synk := make(chan *pixel, 16)
    calcPix := pixGenerator(cam, synk, sem)
    go func() { // generating more goroutines as needed
        for y := 0; y < cam.height; y++ {
            for x := 0; x < cam.width; x++ {
                sem <- struct{}{}
                go calcPix(x, y)
            }
        }
    }()

    img, set := newImg(cam.width, cam.height)
    for i := 0; i < cam.height*cam.width; i++ { // setting pixels of image as they are done
        pix := <- synk
        set(pix.x, pix.y, pix.r, pix.g, pix.b)
    }

    dumpImg(img)
}

func pixGenerator(cam *camera, synk chan *pixel, sem chan struct{}) func(int, int) {
    // fwd := vector(0, 0, -1)
    rht := vector(1, 0, 0)
    up := vector(0, 1, 0)
    topLeft := vector(-float64(cam.width)/2, float64(cam.height)/2, -cam.scrDist)
        
    cam.camPos = vector(0, 0, 0)
    hitSP1 := sphereHitCreate(vector(0, 0, -5), 1) // objects
    groundSP := sphereHitCreate(vector(0, -1005, 0), 1000)    
    
    raysPerPix := 20
    return func(x, y int) { // calculates color of pixel and returns color in a channel
        pix := vector(0, 0, 0)
        for i := 0; i < raysPerPix; i++ {
            randomvec := vector(rand.Float64()-0.5, rand.Float64()-0.5, 0) // we add a random vector to ray to check multiple points within the pixel
            ray := nMatAdd(cam.camPos, topLeft, matScalar(rht, float64(x)), matScalar(up, float64(-y)), randomvec)
            bright1, t1, did1 := hitSP1(cam.camPos, ray)
            bright2, t2, did2 := groundSP(cam.camPos, ray)
            if !did1 && !did2 {
                pix = matAdd(pix, backgroundPix(y, cam.height))
                // set(x, y, round(color[0][0]), round(color[1][0]), round(color[2][0]))
                continue
            }
            var bright float64
            if t1 < t2 && did1 {bright = bright1} else if t2 < t1 && did2 {bright = bright2}
            // if y == 480/2 {fmt.Println(bright, intersectionPoint, t, x, y)}
            // set(x, y, round(absVal(255*bright)), 0, 0)
            pix = matAdd(pix, vector(absVal(255*bright), 0, 0))
        }
        pix = matScalar(pix, 1/float64(raysPerPix)) // submitting pixel
        ixel := pixel{}
        ixel.x, ixel.y = x, y
        ixel.r = int(math.Round((pix[0][0])))
        ixel.g = int(math.Round((pix[1][0])))
        ixel.b = int(math.Round((pix[2][0])))
        synk <- &ixel
        <- sem
    }
}

// returns the color of background. use this if ray dosent hit anything.
// we can use any image instead of just simple stuff.
func backgroundPix(y, height int) [][]float64 {
    t := 1 - float64(y)/float64(height)
    col1 := vector(229, 240, 255)
    col2 := vector(148, 191, 255)
    return matAdd(matScalar(col1, 1-t), matScalar(col2, t))
}

type sphere struct {
    color [][]float64 // color vector
    r float64
    center [][]float64
    material uint8
}

// returns a function that tells if the ray hits this sphere or not
func sphereHitCreate(center [][]float64, r float64) func([][]float64, [][]float64) (float64, float64, bool) {
    return func(camPos, ray [][]float64) (float64, float64, bool) {
        oc := matSub(camPos, center)
        negB := -vecDot(ray, oc)
        bSq := vecDot(ray, ray)
        Dby4 := negB*negB - bSq*(vecDot(oc, oc)-r*r)
        if Dby4 < 0 {return 0, 999999999999, false} // didnt hit
        t := (negB - math.Sqrt(Dby4))/bSq // no +ve sqrt(D) cuz we want min anyway
        if t < 0 {return 0, 999999999999, false} // ray hitting behind the camera
        intersectionPoint := matScalar(ray, t)
        intersectionNormal := matSub(intersectionPoint, center)
        intersectionNormal = vecUnit(intersectionNormal)
        // set(x, y, round(absVal(255*intersectionNormal[0][0])), round(absVal(255*intersectionNormal[1][0])), round(absVal(255*intersectionNormal[2][0])))
        bright := vecDot(intersectionNormal, vecUnit(vector(-1, -1, -2)))
        return bright, t, true
    }
}