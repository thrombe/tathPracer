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
            rayOri := vector(cam.camPos[0][0], cam.camPos[1][0], cam.camPos[2][0])

            shot, depth := shoot(objects, ray, rayOri, 0)
            pix = matAdd(pix, matScalar(shot, 1/float64(depth)))
        }
        submitPix(x, y, matScalar(pix, 1/float64(raysPerPix)), synk)
        <- sem
    }
}

// receieves objects in a list and determines what the given ray hits (ie closest)
func shoot(objects []*sphere, ray, rayOri [][]float64, depth int) ([][]float64, int) {
    depth++
    tmin := 9999999999999999.0
    var color [][]float64
    var hitt *sphere
    for _, ob := range objects {
        bright, t, did := ob.hit(ray, rayOri)
        if !did {continue}
        _ = bright
        // if t < tmin {color, tmin, hitt = matScalar(ob.color, bright), t, ob}
        if t < tmin {color, tmin, hitt = ob.color, t, ob}
    }
    if color == nil {return backgroundPix(ray), depth}
    if depth >= 10 {return color, depth}
    // if rand.Float64() > hitt.albedo {return color, depth} // some light gets absorbed (idk if i should return 0 color or sp.color)
    // if rand.Float64() > hitt.albedo {return vector(0, 0, 0), depth-1}
    // if rand.Float64() > 0.4 {return vector(0, 0, 0), depth}
    ray, rayOri = hitt.getRay(ray, nMatAdd(rayOri, matScalar(ray, 0.999*tmin)))
    shot, dep := shoot(objects, ray, rayOri, depth)
    return matAdd(color, shot), dep
}

// just to improve readability in the other function
func submitPix(x, y int, pix [][]float64, synk chan *pixel) {
    ixel := pixel{}
    ixel.x, ixel.y = x, y
    ixel.r = int(math.Round(255*math.Sqrt(pix[0][0]/255))) // /255 and sqrt *255 for gamma correction
    ixel.g = int(math.Round(255*math.Sqrt(pix[1][0]/255)))
    ixel.b = int(math.Round(255*math.Sqrt(pix[2][0]/255)))
    // ixel.r = int(math.Round(pix[0][0]))
    // ixel.g = int(math.Round(pix[1][0]))
    // ixel.b = int(math.Round(pix[2][0]))
    synk <- &ixel
}

// function that generates some objects. (temporary function) (improves readability in other fn)
func genObjects() []*sphere {
    howMany := 10 // generating random objects
    groundsp := sphere{}
    groundsp.center, groundsp.r, groundsp.color = vector(0, -1005, -5), 1000, vector(0, 255, 0)
    groundsp.material = 0
    objects := make([]*sphere, howMany+4)
    objects[0] = &groundsp
    for i := 1; i < howMany/2+1; i++ { // spheres in front of camera
        sp := sphere{}
        sp.r, sp.color = 1, vector(0, 0, 255)
        sp.material = uint8(math.Round(rand.Float64()))
        // choosing center by randomly distributing spheres in a small area, then displace it up and forward
        // then find the centers wrt goundsp center and multiply the unit vectors by something so they end up on surface of gsp
        // finally add back the gsp.center to displace it back
        vec := matAdd(vector((rand.Float64()-0.5)*30, 0, (rand.Float64()-0.5)*30), vector(0, groundsp.r+sp.r, -40))
        vec = matSub(vec, groundsp.center)
        sp.center = matAdd(matScalar(vec, (groundsp.r+sp.r)/vecSize(vec)), groundsp.center)
        objects[i] = &sp
        // fmt.Println(sp.center)
    }
    for i := howMany/2+1; i < howMany+1; i++ { // spheres behind camera
        sp := sphere{}
        sp.r, sp.color = 1, vector(0, 0, 255)
        sp.material = uint8(math.Round(rand.Float64()))
        vec := matAdd(vector((rand.Float64()-0.5)*30, 0, (rand.Float64()-0.5)*30), vector(0, groundsp.r+sp.r, 10))
        vec = matSub(vec, groundsp.center)
        sp.center = matAdd(matScalar(vec, (groundsp.r+sp.r)/vecSize(vec)), groundsp.center)
        objects[i] = &sp
    }
    sp1 := sphere{}
    sp1.center, sp1.r, sp1.color = vector(0, 0, -5), 1, vector(255, 0, 0)
    sp1.material = 1
    objects[howMany+1] = &sp1
    sp2 := sphere{}
    sp2.center, sp2.r, sp2.color = vector(-2, 0, -5), 1, vector(255, 255, 0)
    sp2.material = 0
    objects[howMany+2] = &sp2
    sp3 := sphere{}
    sp3.center, sp3.r, sp3.color = vector(2, 0, -5), 1, vector(0, 255, 255)
    sp3.material = 2
    objects[howMany+3] = &sp3
    return objects
}

// returns the color of background. use this if ray dosent hit anything.
// we can use any image instead of just simple stuff.
func backgroundPix(ray [][]float64) [][]float64 {
    // t := 1 - float64(y)/float64(height)
    // col1 := vector(229, 240, 255) // whitish
    // col2 := vector(148, 191, 255) // skyish blue
    t := absVal(vecUnit(ray)[1][0])+0.5
    col1 := vector(255, 255, 255) // whitish
    col2 := vector(125, 178, 255) // skyish blue
    return matAdd(matScalar(col1, 1-t), matScalar(col2, t)) // interpolate between both according to y value
}

type sphere struct {
    color [][]float64 // color vector
    r float64
    center [][]float64
    material uint8
    albedo float64 // a degree of how much of the rays hitting the object get absorbed [0, 1)
}

// returns a function that tells if the ray hits this sphere or not
func (sp *sphere) hit(ray, rayOri [][]float64) (float64, float64, bool) {
    oc := matSub(rayOri, sp.center)
    negB := -vecDot(ray, oc)
    bSq := vecDot(ray, ray)
    Dby4 := negB*negB - bSq*(vecDot(oc, oc)-sp.r*sp.r)
    if Dby4 < 0 { // didnt hit
        return  0, 999999999999, false
        }
    t := (negB - math.Sqrt(Dby4))/bSq // no +ve sqrt(D) cuz we want min anyway
    if t < 0 { // ray hitting behind the camera or really close to object ( t < 0.0000001 for really close thing)
        return 0, 999999999999, false
        }
    intersectionPoint := matScalar(ray, t)
    intersectionNormal := matSub(intersectionPoint, sp.center)
    intersectionNormal = vecUnit(intersectionNormal)
    bright := absVal(vecDot(intersectionNormal, vecUnit(vector(-1, -1, -2))))
    return bright, t, true
}

// decides what ray to get. ie reflected or refracted etc based on paameters of the sphere
func (sp *sphere) getRay(ray, point [][]float64) ([][]float64, [][]float64) {
    if sp.material == 0 {
        return sp.diffuse(point)
    } else if sp.material == 1 {
        return sp.reflection(ray, point)
    } else if sp.material == 2 {
        return sp.refraction(ray, point)
    } else {
        return sp.refraction(ray, point)
    }
}

// returns a point on the unit sphere with center 1 unit from the intersection point in the direction of the normal
func (sp *sphere) diffuse(point [][]float64) ([][]float64, [][]float64) {
    unitnormal := vecUnit(matSub(point, sp.center))
    // point = matAdd(point, matScalar(unitnormal, 0.000000001)) // shadow acne ??
    return matAdd(unitnormal, vecUnit(vector(rand.Float64(), rand.Float64(), rand.Float64()))), point
}

// point is the intersection point
func (sp *sphere) reflection(ray, point [][]float64) ([][]float64, [][]float64) {
    normal := matSub(point, sp.center)
    // ray = matScalar(ray, -1)
    rayOri := point
    ray = matAdd(ray, matScalar(normal, -vecDot(ray, vecUnit(normal))*2))
    // ray = matAdd(ray, matScalar(vector(rand.Float64(), rand.Float64(), rand.Float64()), 3)) // fuzzy reflections? didnt work
    return ray, rayOri
}

func (sp *sphere) refraction(ray, point [][]float64) ([][]float64, [][]float64) {
    muglass := 1.3
    muair := 1.0
    radiusvec := matSub(point, sp.center)
    if vecSize(radiusvec) < sp.r {muair, muglass = muglass, muair}
    unitnormal := vecUnit(radiusvec)
    cross := vecCross(ray, unitnormal)
    perplen := -vecSize(cross)*muair/muglass // -ve cuz reasons
    perpdir := vecUnit(vecCross(cross, unitnormal))
    alonglen := -vecDot(ray, unitnormal)
    ray = matAdd(matScalar(unitnormal, -alonglen), matScalar(perpdir, perplen))
    if vecSize(radiusvec) > sp.r {return sp.refraction(ray, matAdd(matScalar(radiusvec, 0.9999), sp.center))}
    return ray, matAdd(matScalar(radiusvec, 1.0001), sp.center)
}

//returns the cross product vectors (i, j, k)
func vecCross(vec1, vec2 [][]float64) [][]float64 {
    return vector(
        vec1[1][0]*vec2[2][0]-vec2[1][0]*vec1[2][0],
        -vec1[0][0]*vec2[2][0]+vec2[0][0]*vec1[2][0],
        vec1[0][0]*vec2[1][0]-vec2[0][0]*vec1[1][0],
    )
}