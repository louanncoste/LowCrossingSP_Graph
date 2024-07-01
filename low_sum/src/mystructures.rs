
use rand::Rng;




pub type Node=u32;

pub struct Edge{
    pub u: Node,
    pub v: Node,
}





//Graph in compressed sparse row format
//Nodes must be numbered with integers from 0 to n-1 
//The Neighborhoods are closed
pub struct Graph{
    deg_sum:Vec<u32>,                //deg_sum[i] contains the number of non-zero element from the first to the i-th line of the adjacency matrix
    adj:Vec<Node>,                   //adj[j] contains the column number of the j-th non-zero element in the adjacency matrix
}

impl Graph{
    pub fn n(&self) -> u32{       //number of vertices
        (self.deg_sum.len()-1) as u32
    }

    pub fn m(&self, undir :bool) -> u32{       //number of edges
        if undir {*self.deg_sum.last().unwrap()/2}    
        else {*self.deg_sum.last().unwrap()}
    }

    pub fn deg(&self, v: Node) -> u32{       
        let begin = self.deg_sum[v as usize];    
        let end=self.deg_sum[(v+1) as usize];
        end-begin
    }

    pub fn neighbors(&self, v: Node) -> &[Node]{
        let begin = self.deg_sum[v as usize] as usize;     
        let end=self.deg_sum[(v+1) as usize] as usize;
        &self.adj[begin..end] 
    }


    //build the graph from a list of edges  
    pub fn from_edges(edges: &Vec<Edge>, undir:bool)-> Self{
        let mut n : Node =0;
        for e in edges{
            if e.u>n {n=e.u;}
            if e.v>n {n=e.v;}
        }
        //add loops
        let mut loops=Vec::with_capacity((n+1) as usize);
        for u in 0..=n{
            loops.push(Edge{u,v:u});
        }

        let mut edges_and_loops =Vec::with_capacity(edges.len()+loops.len());
        for e in edges.iter().chain(loops.iter()){
            edges_and_loops.push(Edge{u:e.u, v:e.v});
        }

        //create deg_sum in 2 parts
        let mut deg_sum: Vec<u32>=vec![0; (n+2)as usize];
        for e in &edges_and_loops{
            deg_sum[e.u as usize]+=1;
            if undir && e.u !=e.v {deg_sum[e.v as usize]+=1;}
        }
        for u in 0..=n{
            deg_sum[(u+1) as usize] +=deg_sum[u as usize]
        }
        //create adj
        let mut adj: Vec<Node> = vec![0; deg_sum[(n+1) as usize] as usize];
        let mut add_adj_fn= |u:Node, v:Node|{
            deg_sum[u as usize]-=1;
            adj[deg_sum[u as usize] as usize]=v;
        };
        for e in edges_and_loops.iter().rev(){
            add_adj_fn(e.u, e.v);
            if undir && e.u!= e.v {add_adj_fn(e.v,e.u);}
        }
        Self{deg_sum, adj}
    }
}





pub struct NodeInfo{
    pub state:u8,                //=0 if not in any path, =1 if one endpoint of a path, =2 if internal in a path (ie. has 2 neighbors in the path)
    pub adjacent_pairs:Vec<u32>, //vector of all pairs that are active in the current iteration. (A pair is represented by its id: its index in active_edge of sampled_pairs)
    pub endpoint:Node,           //If the node is an endpoint of a path (status=1), endpoint represents the node that is the other endpoint of our path 
    pub neigh1:Node,             //Its neighbor when it was first added in the path
    pub neigh2:Node,             //Its other neighbor in the path (the second added)
}


impl NodeInfo{
    pub fn create_node_info(i:u32)->Self{
        let adjacent_pairs=Vec::new();
        Self{state:0, adjacent_pairs, endpoint:i, neigh1:i, neigh2:i}
    }

    pub fn init_adjacent_pairs(&mut self){
        self.adjacent_pairs=Vec::new();
    }
}





pub struct Pair{
    pub u:Node,
    pub v:Node,
    pub id:u32,
    pub index_in_cost:u32,
}



pub struct MyPairs{
    pub active_pairs:Vec<Pair>,     //contains all the pairs we sampled at the beginning of compute_path_iteration
    pub cost_vect:Vec<u32>,         //sorted by decreasing order
}


impl MyPairs{

    //Create a MyPairs element, sampling 'sample_size_per_vertex' edges per vertex
    pub fn create_sample_pairs(nodelist: &Vec<Node>, sample_size_per_vertex:usize, all_node_info: &mut Vec<NodeInfo>)->Self{
        let mut active_pairs:Vec<Pair>=Vec::with_capacity(nodelist.len()*sample_size_per_vertex);
        let mut cp=0;
        for u in nodelist{
            for _ in 0..sample_size_per_vertex{
                let v= loop{
                    let k=rand::thread_rng().gen_range(0..nodelist.len());
                    if nodelist[k] != *u && nodelist[k] != all_node_info[*u as usize].endpoint{
                        break nodelist[k];
                    }
                };
                active_pairs.push(Pair{u:*u, v, id:cp, index_in_cost:0});
                all_node_info[*u as usize].adjacent_pairs.push(cp);
                all_node_info[v as usize].adjacent_pairs.push(cp);
                cp+=1;
            }
        }
        Self{active_pairs, cost_vect: Vec::new()}
    }


    /* 
    pub fn create_sample_pairs_from_edgelist(edgelist: &Vec<Edge>, all_node_info: &mut Vec<NodeInfo>)->Self{
        let mut active_pairs:Vec<Pair>=Vec::with_capacity(edgelist.len());
        let mut weight_buckets:Vec<Vec<u32>>=vec![Vec::with_capacity(edgelist.len())];
        let is_neigh_x:Vec<bool>=vec![false; all_node_info.len()];
        let mut cp=0;
        for e in edgelist{
            if e.u == e.v {continue;}
            active_pairs.push(Pair{u:e.u, v:e.v, id:cp, intersection_count:0, index_in_my_weight_bucket:cp});
            weight_buckets[0].push(cp);
            all_node_info[e.u as usize].adjacent_pairs.push(cp);
            all_node_info[e.v as usize].adjacent_pairs.push(cp);
            cp+=1;
        }
        Self{active_pairs, weight_buckets, min_bucket_index:0, is_neigh_x, capacity:cp}
    }*/


   // pub fn resample_pairs()

    pub fn compute_costs(& mut self, g:&Graph){  
        let mut all_crossing_numbers: Vec<u32>=vec![0; self.active_pairs.len()];
        let mut sorted_by_cost=Vec::with_capacity(self.active_pairs.len());
       
        for p in &self.active_pairs{
            sorted_by_cost.push(p.id);
            let u=p.u;
            let v=p.v;
            let mut sum=0; //sym difference of neighborhoods of u and v (=crossing number of p)

            let neigh_u = g.neighbors(u); //those are sorted
            let neigh_v = g.neighbors(v);

            let mut i: usize = 0;
            let mut j: usize = 0;
            while i < neigh_u.len() && j < neigh_v.len() {
                if neigh_u[i] < neigh_v[j] {
                    sum += 1;
                    i += 1;
                } else if neigh_u[i] > neigh_v[j] {
                    sum += 1;
                    j += 1;
                } else {
                    i += 1;
                    j += 1;
                }
            }
            sum += neigh_u.len() - i;
            sum += neigh_v.len() - j;
            
            all_crossing_numbers[p.id as usize]=sum as u32;
        } 

        sorted_by_cost.sort_by(|&a, &b| all_crossing_numbers[b as usize].cmp(&all_crossing_numbers[a as usize])); //sorted in decreasing order
        for ind in 0..sorted_by_cost.len(){
            self.active_pairs[sorted_by_cost[ind] as usize].index_in_cost=ind as u32;
        }
        self.cost_vect=sorted_by_cost;
    }


    pub fn get_minimal_cost_pair(&mut self, all_node_info: &Vec<NodeInfo>)-> &Pair{
        loop{
            if self.cost_vect.len()==0 {return &Pair{u:u32::MAX, v:0, id:0, index_in_cost:0}}
            let p_indx= *self.cost_vect.last().unwrap();
            let u= self.active_pairs[p_indx as usize].u;
            let v= self.active_pairs[p_indx as usize].v;
            if all_node_info[u as usize].endpoint == v || all_node_info[u as usize].state >=2 || all_node_info[v as usize].state >=2{
                self.cost_vect.pop();
                self.active_pairs[p_indx as usize].index_in_cost = u32::MAX;
            }
            else{
                break &self.active_pairs[p_indx as usize];
            }
        }
    }




}


