use std::cmp;
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
    pub intersection_count:u32,
    pub index_in_my_weight_bucket:u32,
}



pub struct MyPairs{
    pub active_pairs:Vec<Pair>,     //contains all the pairs we sampled at the beginning of compute_path_iteration
    weight_buckets:Vec<Vec<u32>>,  //vector of 'buckets'. The 'bucket' at weight_buckets[i] contains the indices in active_pairs of all the pairs of weight 1/(2**i)
    min_bucket_index:u32,           //the index of the first non-empty bucket
    is_neigh_x: Vec<bool>,          //is useful in update_weights_pairs
}


impl MyPairs{

    //Create a MyPairs element, sampling 'sample_size_per_vertex' edges per vertex
    pub fn create_sample_pairs(nodelist: &Vec<Node>, sample_size_per_vertex:usize, all_node_info: &mut Vec<NodeInfo>)->Self{
        let mut active_pairs:Vec<Pair>=Vec::with_capacity(nodelist.len()*sample_size_per_vertex);
        let mut weight_buckets:Vec<Vec<u32>>=vec![Vec::with_capacity(nodelist.len()*sample_size_per_vertex)];
        let is_neigh_x:Vec<bool>=vec![false; all_node_info.len()];
        let mut cp=0;
        for u in nodelist{
            for _ in 0..sample_size_per_vertex{
                let v= loop{
                    let k=rand::thread_rng().gen_range(0..nodelist.len());
                    if nodelist[k] != *u && nodelist[k] != all_node_info[*u as usize].endpoint{
                        break nodelist[k];
                    }
                };
                active_pairs.push(Pair{u:*u, v, id:cp, intersection_count:0, index_in_my_weight_bucket:cp});
                weight_buckets[0].push(cp);
                all_node_info[*u as usize].adjacent_pairs.push(cp);
                all_node_info[v as usize].adjacent_pairs.push(cp);
                cp+=1;
            }
        }
        Self{active_pairs, weight_buckets, min_bucket_index:0, is_neigh_x}
    }

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
        Self{active_pairs, weight_buckets, min_bucket_index:0, is_neigh_x}
    }


    // samples one pair according to current weights
    pub fn weighted_sample_pair(&mut self, all_node_info: &Vec<NodeInfo>,pair_sampling_search_radius:usize)->&Pair{
        if self.weight_buckets.len()==0 || self.min_bucket_index >= self.weight_buckets.len() as u32{
            return &Pair{u: u32::MAX, v:0, id:0, intersection_count:0, index_in_my_weight_bucket:0}
        }
        loop{    
            //until find a pair e that is not two endpoints of a same path
            let e_indx= 'sampling: loop{
                //until find a pair e in our search radius
                let max_pair_index= cmp::min(self.min_bucket_index as usize+pair_sampling_search_radius, self.weight_buckets.len()-1);
                let mut total_weight=0;
                let mut weight=1;
                for i in (self.min_bucket_index as usize..=max_pair_index).rev(){
                    total_weight+=self.weight_buckets[i].len()*weight;
                    weight<<=1;                //weights are multiple of 2
                    if i==0 {break};
                }

                weight>>=1;

                if total_weight==0 {
                    return &Pair{u: u32::MAX, v:0, id:0, intersection_count:0, index_in_my_weight_bucket:0}
                }
                let r=rand::thread_rng().gen_range(0..total_weight);
                let mut curr_sum=0;
                for i in self.min_bucket_index as usize..=max_pair_index{
                    curr_sum+=self.weight_buckets[i].len()*weight;
                    if curr_sum>r {
                        let rd =rand::thread_rng().gen_range(0..self.weight_buckets[i].len());
                        break 'sampling  self.weight_buckets[i][rd] as usize;
                    }
                    weight>>=1;
                }
            };
            let u:Node = self.active_pairs[e_indx].u;
            let v:Node = self.active_pairs[e_indx].v;
            if all_node_info[u as usize].endpoint==v{
                self.delete_pair(e_indx as u32);
            }
            else{
                break &self.active_pairs[e_indx];  
            }     
        }
    }


    //a usual swap_remove for the vector weight_buckets (but also handles the index field)
    fn swap_remove_weight_bucket(&mut self, weight:usize, index_bucket:usize){
        if index_bucket== self.weight_buckets[weight].len()-1 {
            self.weight_buckets[weight].pop();
        }
        else{
            self.weight_buckets[weight].swap_remove(index_bucket);
            self.active_pairs[self.weight_buckets[weight][index_bucket]as usize].index_in_my_weight_bucket=index_bucket as u32;
        }
    } 


    //increment the weight of the pair 'pair_id' in the weight_bucket
    fn increment_weight_p(&mut self, pair_id: u32){

        let mypair_intersect= self.active_pairs[pair_id as usize].intersection_count as usize;
        let mypair_ind_weightbucket= self.active_pairs[pair_id as usize].index_in_my_weight_bucket as usize;       

        self.swap_remove_weight_bucket(mypair_intersect ,mypair_ind_weightbucket);

        //If we were in the max bucket, we create a new one
        if self.weight_buckets.len() == mypair_intersect+1{
            self.weight_buckets.push(Vec::with_capacity(self.active_pairs.len()));
        }

        self.weight_buckets[mypair_intersect+1].push(pair_id);

        //update the indices
        let pair=&mut self.active_pairs[pair_id as usize];
        pair.index_in_my_weight_bucket=(self.weight_buckets[mypair_intersect+1].len()-1) as u32;
        pair.intersection_count+=1;


        //if we were the last one of the min bucket, we deallocate the memory
        if self.weight_buckets[self.min_bucket_index as usize].is_empty(){
            self.weight_buckets[self.min_bucket_index as usize]=Vec::new();
            self.min_bucket_index+=1;
        }
    }


    //delete a pair from weight_buckets, but not from active_pairs
    pub fn delete_pair(&mut self, pair_id : u32){
        let mypair_intersect= self.active_pairs[pair_id as usize].intersection_count ;
        let mypair_ind_weightbucket= self.active_pairs[pair_id as usize].index_in_my_weight_bucket ;  
        
        if mypair_ind_weightbucket==u32::MAX {return;}

        self.swap_remove_weight_bucket(mypair_intersect as usize,mypair_ind_weightbucket as usize);


        let pair = &mut self.active_pairs[pair_id as usize];
        pair.index_in_my_weight_bucket=u32::MAX;
        
        while (self.min_bucket_index<self.weight_buckets.len() as u32) && (self.weight_buckets[self.min_bucket_index as usize].is_empty()){
            self.min_bucket_index+=1;
        }
    }


    //increment weight of all 'active' pairs crossing x
    pub fn update_weights_pairs(&mut self, x : &Neighborhood, g: &Graph, all_node_info: &mut Vec<NodeInfo>) {
        //for all neighbors of x (itself included), we check if all the active pairs from those neighbors are crossing or not
        //If they do, we increment their weight 

        self.is_neigh_x.fill(false);
        let all_x_neigh=g.neighbors(x.id);
        for node in all_x_neigh{
            self.is_neigh_x[*node as usize]=true;
        }

        for node in all_x_neigh{
            let mut i=0; 
            let mut adj_len=all_node_info[*node as usize].adjacent_pairs.len();
            while i< adj_len{                     //if node is not in nodelist, its adjacent_pairs is empty
                let p_indx=all_node_info[*node as usize].adjacent_pairs[i];
                let p= &mut self.active_pairs[p_indx  as usize];
                if p.index_in_my_weight_bucket==u32::MAX{
                    all_node_info[*node as usize].adjacent_pairs.swap_remove(i);
                    adj_len-=1;
                }
                else{
                    i+=1;
                    if all_node_info[p.u as usize].endpoint !=p.v{
                        let other_node= if p.u== *node {p.v} else {p.u};
                        if !self.is_neigh_x[other_node as usize]{
                            self.increment_weight_p(p_indx as u32);
                        }
                    }
                    else {   //u and v are endpoint of a same path
                        self.delete_pair(p_indx);
                    }
                }   
            }
        }
    }

}




pub struct Neighborhood{
    pub id:Node,
    pub intersection_count:u32,
    pub index_in_my_weight_bucket:u32,
}

pub struct MyNeighborhoods{
    all_neighborhoods:Vec<Neighborhood>,  //all_neighborhoods[i] contains the Neighborhood of id 'i'. All vertices of the initial graph are contained in this vector
    weight_buckets:Vec<Vec<u32>>,  //vector of 'buckets'. The 'bucket' at weight_buckets[i] contains the indices in 'all_neighborhoods' of all the neighborhoods of weight 2**i
    is_crossing: Vec<u8>,   //useful in update_weights_neighborhoods

}



impl MyNeighborhoods{

    pub fn create_neighborhoods(n:usize)->Self{
        let mut all_neighborhoods: Vec<Neighborhood>=Vec::with_capacity(n);
        let mut weight_buckets:Vec<Vec<u32>>=vec![Vec::with_capacity(n)];
        let is_crossing:Vec<u8>=vec![0; n];
        for i in 0..n as u32{
            all_neighborhoods.push(Neighborhood{id:i, intersection_count:0, index_in_my_weight_bucket:i});
            weight_buckets[0].push(i);
        }
        Self{all_neighborhoods, weight_buckets, is_crossing}
    }


    //init weight_buckets with only the nodes of nodelist
    pub fn init_weights(&mut self, nodelist :&Vec<Node>){
        let mut new_weight_buckets :Vec<Vec<u32>>=vec![Vec::with_capacity(nodelist.len())];
        for i in 0..nodelist.len(){
            let n_indx=nodelist[i] ;
            new_weight_buckets[0].push(n_indx);
            self.all_neighborhoods[n_indx as usize].intersection_count=0;
            self.all_neighborhoods[n_indx as usize].index_in_my_weight_bucket=i as u32;
        }
        self.weight_buckets=new_weight_buckets;
    }


    // samples one neighborhood according to current weights
    pub fn weighted_sample_neighborhood(&self, neighborhood_sampling_search_radius:usize)->&Neighborhood{
        'sampling: loop{
            let mut total_weight=0;
            let mut weight=1;
            let max_bucket_index=self.weight_buckets.len()-1;
            let min_neigh_index =if neighborhood_sampling_search_radius< max_bucket_index {max_bucket_index-neighborhood_sampling_search_radius} else{0};
            for i in min_neigh_index..=max_bucket_index{
                total_weight+=self.weight_buckets[i].len()*weight;
                weight<<=1;
            }
            weight>>=1;
            let r=rand::thread_rng().gen_range(0..total_weight);
            let mut curr_sum=0;
            for i in (min_neigh_index..=max_bucket_index).rev(){
                curr_sum+=self.weight_buckets[i].len()*weight;
                if curr_sum>r {
                    let rd =rand::thread_rng().gen_range(0..self.weight_buckets[i].len());
                    break 'sampling  &self.all_neighborhoods[self.weight_buckets[i][rd] as usize];
                }
                weight>>=1;
                if i==0 {break;}
            }
        }
    }


    fn swap_remove_weight_bucket(&mut self, weight:usize, index_bucket:usize){
        if index_bucket== self.weight_buckets[weight].len()-1 {
            self.weight_buckets[weight].pop();
        }
        else{
            self.weight_buckets[weight].swap_remove(index_bucket);
            self.all_neighborhoods[self.weight_buckets[weight][index_bucket] as usize].index_in_my_weight_bucket=index_bucket as u32;
        }
    } 


    //increment the weight of the neighborhood 'neigh_id' in the weight_bucket
    fn increment_weight_n(&mut self, neigh_id: &u32){

        let myneigh_intersect= self.all_neighborhoods[*neigh_id as usize].intersection_count as usize;
        let myneigh_ind_weightbucket= self.all_neighborhoods[*neigh_id as usize].index_in_my_weight_bucket as usize;       

        self.swap_remove_weight_bucket(myneigh_intersect,myneigh_ind_weightbucket);

        //If we were in the max bucket, we create a new one
        if self.weight_buckets.len() == myneigh_intersect+1{
            self.weight_buckets.push(Vec::with_capacity(self.all_neighborhoods.len()));
        }

        self.weight_buckets[myneigh_intersect+1].push(*neigh_id);

        //update the indices
        let neigh=&mut self.all_neighborhoods[*neigh_id as usize];
        neigh.index_in_my_weight_bucket=(self.weight_buckets[myneigh_intersect+1].len()-1) as u32;
        neigh.intersection_count+=1;
    }
    


    //we only delete from weight_buckets
    pub fn delete_neighborhood(&mut self, neigh_id : u32){
        let myneigh_intersect= self.all_neighborhoods[neigh_id as usize].intersection_count;
        let myneigh_ind_weightbucket= self.all_neighborhoods[neigh_id as usize].index_in_my_weight_bucket;   

        if myneigh_ind_weightbucket==u32::MAX {return;}

        self.swap_remove_weight_bucket(myneigh_intersect as usize,myneigh_ind_weightbucket as usize);

        let neigh=&mut self.all_neighborhoods[neigh_id as usize];
        neigh.index_in_my_weight_bucket=u32::MAX;

    }
    

    //increment weight of all the neighborhoods crossed by the pair 'e' 
    pub fn update_weights_neighborhoods(&mut self, e : &Pair, g: &Graph){
        self.is_crossing.fill(0); 
        //is_crossing[i] = 0 (e is outside neighborhood i), 1 (e crosses neighborhood i), or 2 (e is inside neighborhood i)

        for node in [e.u, e.v]{
            for node_neigh in g.neighbors(node){
                self.is_crossing[*node_neigh as usize]+=1;
            }
        }

        for node in [e.u, e.v]{
            for node_neigh in g.neighbors(node){
                if self.all_neighborhoods[*node_neigh as usize].index_in_my_weight_bucket==u32::MAX {continue;}
                if self.is_crossing[*node_neigh as usize]==1{
                    self.increment_weight_n(node_neigh);
                }
            }
        }
    }

}