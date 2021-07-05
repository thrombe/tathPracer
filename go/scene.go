package main
import (
	"math"
	"math/rand"
)

// function that generates some objects. (temporary function) (improves readability in other fn)
func genObjects() []*sphere {
    howMany := 10 // generating random objects
    groundsp := sphere{}
    groundsp.center, groundsp.r, groundsp.color = vector(0, -1005, -5), 1000, vector(0, 255, 0)
    groundsp.material = 0
    objects := make([]*sphere, howMany+5)
    objects[0] = &groundsp
    for i := 1; i < howMany/2+1; i++ { // spheres in front of camera
        sp := sphere{}
        sp.r, sp.color = 1, vector(0, 0, 255)
        sp.material = uint8(math.Round(rand.Float64()*2.5))
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
    sp1.material = 2
    objects[howMany+1] = &sp1
    sp2 := sphere{}
    sp2.center, sp2.r, sp2.color = vector(-2, 0, -5), 1, vector(255, 255, 0)
    sp2.material = 0
    objects[howMany+2] = &sp2
    sp4 := sphere{}
    sp4.center, sp4.r, sp4.color = vector(1, 0.7, -30), 1, vector(255, 0, 255)
    sp4.material = 0
    objects[howMany+4] = &sp4
    sp3 := sphere{}
    sp3.center, sp3.r, sp3.color = vector(2, 0, -5), 1, vector(0, 255, 255)
    sp3.material = 1
    objects[howMany+3] = &sp3
    return objects
}