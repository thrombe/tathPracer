package main
import (
	"math"
)

type pixel struct {
    x, y int
    r, g, b int
}

// main func that calls other funcs and handles concurrency
func engine() {
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
