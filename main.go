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

type pixel struct {
    x, y int
    r, g, b int
}

// main func that calls other funcs and handles concurrency
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

// returns a func that shoots multiple rays per pixel to figure out the color of that pixel
func pixGenerator(cam *camera, synk chan *pixel, sem chan struct{}) func(int, int) {
    // fwd := vector(0, 0, -1)
    rht := vector(1, 0, 0)
    up := vector(0, 1, 0)
    topLeft := vector(-float64(cam.width)/2, float64(cam.height)/2, -cam.scrDist)
    cam.camPos = vector(0, 0, 0)

    objects := genObjects()
    
    raysPerPix := 2 // 20
    return func(x, y int) { // calculates color of pixel and returns color in a channel
        pix := vector(0, 0, 0)
        for i := 0; i < raysPerPix; i++ {
            randomvec := vector(rand.Float64()-0.5, rand.Float64()-0.5, 0) // we add a random vector to ray to check multiple points within the pixel
            ray := nMatAdd(cam.camPos, topLeft, matScalar(rht, float64(x)), matScalar(up, float64(-y)), randomvec)

            pix = matAdd(pix, shoot(objects, ray, cam, y))
        }
        submitPix(x, y, matScalar(pix, 1/float64(raysPerPix)), synk)
        <- sem
    }
}

// receieves objects in a list and determines what the given ray hits (ie closest)
func shoot(objects []*sphere, ray [][]float64, cam *camera, y int) [][]float64 {
    tmin := 9999999999999999.0
    var color [][]float64
    for _, ob := range objects {
        bright, t, did := ob.hit(cam, ray)
        if !did {continue}
        if t < tmin {color, tmin = matScalar(ob.color, bright), t}
    }
    if color == nil {color = backgroundPix(y, cam.height)}
    return color
}

// just to improve readability in the other function
func submitPix(x, y int, pix [][]float64, synk chan *pixel) {
    ixel := pixel{}
    ixel.x, ixel.y = x, y
    ixel.r = int(math.Round((pix[0][0])))
    ixel.g = int(math.Round((pix[1][0])))
    ixel.b = int(math.Round((pix[2][0])))
    synk <- &ixel
}

// function that generates some objects. (temporary function)
func genObjects() []*sphere {
    howMany := 5 // generating random objects
    groundsp := sphere{}
    groundsp.center, groundsp.r, groundsp.color = vector(0, -1005, -5), 1000, vector(0, 255, 0)
    objects := make([]*sphere, howMany+2)
    objects[0] = &groundsp
    for i := 1; i < howMany+1; i++ {
        sp := sphere{}
        sp.r, sp.color = 1, vector(0, 0, 255)
        // choosing center by randomly distributing spheres in a small area, then displace it up and forward
        // then find the centers wrt goundsp center and multiply the unit vectors by something so they end up on surface of gsp
        // finally add back the gsp.center to displace it back
        vec := matAdd(vector((rand.Float64()-0.5)*30, 0, (rand.Float64()-0.5)*30), vector(0, groundsp.r+sp.r, -40))
        vec = matSub(vec, groundsp.center)
        sp.center = matAdd(matScalar(vec, (groundsp.r+sp.r)/vecSize(vec)), groundsp.center)
        objects[i] = &sp
        // fmt.Println(sp.center)
    }
    sp1 := sphere{}
    sp1.center, sp1.r, sp1.color = vector(0, 0, -5), 1, vector(255, 0, 0)
    objects[howMany+1] = &sp1
    return objects
}

// returns the color of background. use this if ray dosent hit anything.
// we can use any image instead of just simple stuff.
func backgroundPix(y, height int) [][]float64 {
    t := 1 - float64(y)/float64(height)
    col1 := vector(229, 240, 255) // whitish
    col2 := vector(148, 191, 255) // skyish blue
    return matAdd(matScalar(col1, 1-t), matScalar(col2, t)) // interpolate between both according to y value
}

type sphere struct {
    color [][]float64 // color vector
    r float64
    center [][]float64
    material uint8
}

// returns a function that tells if the ray hits this sphere or not
func (sp *sphere) hit(cam *camera, ray [][]float64) (float64, float64, bool) {
    oc := matSub(cam.camPos, sp.center)
    negB := -vecDot(ray, oc)
    bSq := vecDot(ray, ray)
    Dby4 := negB*negB - bSq*(vecDot(oc, oc)-sp.r*sp.r)
    if Dby4 < 0 { // didnt hit
        return  0, 999999999999, false
        }
    t := (negB - math.Sqrt(Dby4))/bSq // no +ve sqrt(D) cuz we want min anyway
    if t < 0 { // ray hitting behind the camera
        return 0, 999999999999, false
        }
    intersectionPoint := matScalar(ray, t)
    intersectionNormal := matSub(intersectionPoint, sp.center)
    intersectionNormal = vecUnit(intersectionNormal)
    bright := absVal(vecDot(intersectionNormal, vecUnit(vector(-1, -1, -2))))
    return bright, t, true
}