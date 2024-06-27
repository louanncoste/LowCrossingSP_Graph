//get the created Spanning Path and all the crossing numbers
fn get_crossing_nb_sum(spanning_path_order: &Vec<u32>, g:&Graph){  //similar to the function update_weights_neighborhoods
    let mut all_crossing_numbers: Vec<u32>=vec![0; g.n()];
    
    let mut i=0;
    while i< g.n()-1{
        let u=spanning_path_order[i];
        let v=spanning_path_order[i+1];

        let mut is_crossing =vec![0;g.n()]; 
        for node in [u, v]{
            for node_neigh in g.neighbors(node){
                is_crossing[*node_neigh as usize]+=1;
            }
        }

        for node in [u, v]{
            for node_neigh in g.neighbors(node){
                if is_crossing[*node_neigh as usize]==1{
                    all_crossing_numbers[*node_neigh as usize]+=1;
                }
            }
        }
        i+=1;
    }
    all_crossing_numbers.iter().sum()
}


/* 

pub type Distance = u32;

pub struct BFS {
    pub dist: Vec<Distance>,
    pub fifo: Vec<Node>,
    front: usize, // index of the front of the fifo
}

impl BFS {

    const INFINITY: Distance = Distance::MAX;

    pub fn new(n: usize) -> Self {
        let dist = vec![Self::INFINITY; n];
        let fifo = Vec::with_capacity(n);
        BFS { dist, fifo, front: 0 }
    }

    pub fn clear(&mut self) {
        self.dist.fill(Self::INFINITY);
        self.fifo.clear();
        self.front = 0;
    }

    pub fn dist(&self, v: Node) -> Distance {
        self.dist[v as usize]
    }

    // return the number of nodes visited so far
    pub fn bfs<G: AGraph>(&mut self, g: &G, v: Node) -> usize {
        self.check_size(g.n());
        self.dist[v as usize] = 0;
        self.fifo.push(v);
        while self.front < self.fifo.len() {
            let v = self.fifo[self.front];
            //eprintln!("bfs: {v}");
            self.front += 1;
            let d = self.dist[v as usize] + 1;
            for &w in g.neighbors(v) {
                if self.dist[w as usize] == Self::INFINITY {
                    self.dist[w as usize] = d;
                    self.fifo.push(w);
                }
            }
        }
        self.front
    }

    fn check_size(&self, n: usize) {
        if n > self.dist.len() { 
            panic!("Graph of {} nodes is too large for this BFS struct of size {}!", 
                   n, self.dist.len()) 
        }
    }

    // Connected Components as a vector indicating beginning index in fifo.
    // The length of the vector is the number of connected components + 1.
    pub fn cc<G: AGraph>(&mut self, g: &G) -> Vec<usize> {
        self.check_size(g.n());
        let mut cc: Vec<usize> = Vec::with_capacity(g.n()+1);
        cc.push(0);
        for v in 0..g.n() as Node {
            if self.dist[v as usize] == Self::INFINITY {
                let iend = self.bfs(g, v); // end of cc of v
                cc.push(iend);
            }
        }
        cc
    }

    pub fn largest_cc(&self, cc : &Vec<usize>) -> &[Node] {
        let mut c_max = 0;
        let mut sz_max = cc[c_max+1];
        for c in 1 .. cc.len() - 1 {
            let sz = cc[c+1] - cc[c];
            if sz > sz_max {
                sz_max = sz;
                c_max = c;
            }
        }
        let b = cc[c_max];
        let e = cc[c_max+1];
        &self.fifo[b..e]
    }

}




fn bfs_numbering(&self) -> Self {
    let mut bfs = BFS::new(self.n());
    bfs.cc(self);
    assert_eq!(bfs.fifo.len(), self.n());
    let mut perm = vec![0 as Node; self.n()];
    for (i, &v) in bfs.fifo.iter().enumerate() {
        perm[v as usize] = i as Node;
    }
    self.numbering(&perm)
}
*/

























fn main(){
    let args: Vec<String> = env::args().collect();

    //get the graph
    let graph_path = &args[1];
    let graph_string = fs::read_to_string(file_path)
        .expect("Should have been able to read the file");
    let mut edges: Vec<Edge>= Vec::new();
    for line in graph_string.lines() {
        let mut parts = line.split_whitespace(); // Split by whitespace
        let u: Node = parts.next().expect("Invalid edge format").parse().expect("Invalid number");
        let v: Node = parts.next().expect("Invalid edge format").parse().expect("Invalid number");
        let edg=Edge{u,v};
        edges.push(edg);
    }
    let g=Graph::from_edges(&mut edges, true);

    let given_order=vec![0;g.n()];
    for i in 0..g.n(){
        given_order[i]=i;
    }


    let low_crossing_order=vec![0;g.n()];
    let lowcross_path = &args[2];
    let lowcross_string = fs::read_to_string(file_path)
        .expect("Should have been able to read the file");
    let mut i=0;
    for line in lowcross_string.lines() {
        let u: Node = line.next().expect("Invalid node format").parse().expect("Invalid number");
        low_crossing_order[i]=u;
        i+=1;
    }

    //let bfs_order=get_bfs(&g);

    println!("Given order Crossing Number Sum= {}", get_crossing_nb_sum(&given_order, &g));
    println!("LowCross order Crossing Number Sum= {}", get_crossing_nb_sum(&low_crossing_order, &g));

}