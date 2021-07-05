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

func vecLerp(vec1, vec2 [][]float64, t float64) [][]float64 {
    return matAdd(matScalar(vec1, t), matScalar(vec2, 1-t))
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

            shot := shoot(objects, ray, rayOri, 0, vector(0, 0, 0))
            // pix = matAdd(pix, matScalar(shot, 1/float64(depth)))
            pix = matAdd(pix, shot)
        }
        submitPix(x, y, matScalar(pix, 1/float64(raysPerPix)), synk)
        <- sem
    }
}

// receieves objects in a list and determines what the given ray hits (ie closest)
func shoot(objects []*sphere, ray, rayOri [][]float64, depth int, color [][]float64) ([][]float64) {
    depth++
    tmin := 9999999999999999.0
    var hitt *sphere
    for _, ob := range objects {
        t, did := ob.hit(ray, rayOri)
        if !did {continue}
        if t < tmin {tmin, hitt = t, ob}
    }
    if hitt == nil {return backgroundPix(ray)}
    if depth >= 10 {return vector(0, 0, 0)} // if its bouncing around too much without end, return 0 color ?
    color = hitt.getCol(color)
    // color = hitt.color
    ray, rayOri = hitt.getRay(ray, rayOri, tmin)
    shot := shoot(objects, ray, rayOri, depth, color)
    color = vecLerp(color, shot, 0.5) // do this inside sp.getCol
    return color
}