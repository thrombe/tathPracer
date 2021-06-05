package main
import (
	"math/rand"
)

//returns the cross product vectors (i, j, k)
func vecCross(vec1, vec2 [][]float64) [][]float64 {
    return vector(
        vec1[1][0]*vec2[2][0]-vec2[1][0]*vec1[2][0],
        -vec1[0][0]*vec2[2][0]+vec2[0][0]*vec1[2][0],
        vec1[0][0]*vec2[1][0]-vec2[0][0]*vec1[1][0],
    )
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