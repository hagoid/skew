fetch('data/18.json')
    .then((response) => response.json())
    .then((json) => console.log((<Array<Array<number>>> json).map(x => create(x, 18))));

const helloWorld = ()=> {
    const element = document.getElementById("hello-world");
    if (element) {
        element.textContent = "Hello, World! Ja som hago";
    }
}

let filename = "data/18.json"

type SkewMorphism = {
    series: Series,
    n: number
}

let skews: Array<SkewMorphism> = []

type Series = Array<number>
type Sequence = Array<number>
type Orbit = Array<number>
type Orbits = Array<Orbit>

const create = (series: Series, n: number): SkewMorphism => {
    return {series: series, n: n}
}


const read = (filename: string) => {

}

helloWorld();