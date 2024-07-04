pub mod mystructures;

use std::fs;
use std::env;
use std::cmp;
use std::fs::{File, OpenOptions};
use std::io::{Write, BufWriter};
use std::time::Instant;
use std::str::FromStr;



use mystructures::Node;
use mystructures::Edge;
use mystructures::Graph;
use mystructures::NodeInfo;
use mystructures::Pair;
use mystructures::MyPairs;
use mystructures::Neighborhood;
use mystructures::MyNeighborhoods;




//get the created Spanning Path and all the crossing numbers
fn get_results(all_node_info: &Vec<NodeInfo>, spanning_path: &mut Vec<Node>, all_crossing_numbers: &mut Vec<u32>, begin: Node, end: Node, g:&Graph){
    let mut current_node=begin;
    let mut previous_node=begin;
    let mut next_node;
    for i in 0..all_node_info.len()-1{
        next_node=if all_node_info[current_node as usize].neigh1==previous_node{all_node_info[current_node as usize].neigh2} else{all_node_info[current_node as usize].neigh1};
        spanning_path[i]=current_node;
        update_crossing(current_node, next_node, all_crossing_numbers,g);

        previous_node=current_node;
        current_node=next_node;
    }
    assert!(current_node==end);
    spanning_path[all_node_info.len()-1]=current_node;
}


//aux function of get_results, for an edge (u,v), we update the crossing of all neighborhoods crossed by (u,v)
fn update_crossing(u:Node, v:Node, all_crossing_numbers: &mut Vec<u32>, g:&Graph){  //similar to the function update_weights_neighborhoods
    let mut is_crossing =vec![0;all_crossing_numbers.len()]; 
    //is_crossing[i] = 0 :(u,v) is outside neighborhood i, 1:(u,v) crosses neighborhood i, or 2:(u,v) is inside neighborhood i

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
}








//////////////////SPANNING PATH ALGORITHM///////////////////




//Auxiliary function 
//Called several times by the main function 'compute_path' in order to regularly reset the weights 
fn compute_path_iteration(nodelist: &Vec<Node>, g: &Graph, all_node_info: &mut Vec<NodeInfo>, neighborhoods: &mut MyNeighborhoods, alpha: f32, it_nb:u32, edgelist: &Vec<Edge>){
   
    neighborhoods.init_weights(&nodelist);
    for node in &mut *all_node_info {
        node.init_adjacent_pairs();
    }

    let sample_size_per_vertex= (nodelist.len() as f32).powf(alpha) as usize +1;

    //println!("edgelist len= {}, sample_size_per_vertex*n={}", edgelist.len(), sample_size_per_vertex*nodelist.len());
    let mut sampled_pairs:MyPairs;
    if it_nb>0 {
        sampled_pairs=MyPairs::create_sample_pairs(nodelist, sample_size_per_vertex, all_node_info);
    }

    else {
        sampled_pairs=MyPairs::create_sample_pairs_from_edgelist(edgelist, all_node_info);
    }

    let nb_real_iteration: u32 = cmp::max(1,(nodelist.len()as f32* 0.5)as u32);
    let nb_warmup= nb_real_iteration/2 as u32; //XXTODO


    for i in 0.. nb_warmup+nb_real_iteration{     
        if i==nb_warmup {   // initalize plane weights after warm-up
            neighborhoods.init_weights(&nodelist);
        }
        let log2n= f32::log2(nodelist.len()as f32);
        let pair_sampling_search_radius= (2.0* log2n +1.0) as usize;
        let neighborhood_sampling_search_radius=(2.0*log2n +1.0) as usize; 
        
        let x: &Neighborhood = neighborhoods.weighted_sample_neighborhood(neighborhood_sampling_search_radius); 
        sampled_pairs.update_weights_pairs(&x, g, all_node_info);

        let e: &Pair= sampled_pairs.weighted_sample_pair(all_node_info, pair_sampling_search_radius); 
        let u=e.u;
        let v=e.v;
        if u==u32::MAX {return;} // if sampled_pairs is empty 
        neighborhoods.update_weights_neighborhoods(&e, g);


        //After warm-up, we start adding the pair e to our spanning path
        if i>= nb_warmup {

            for (node,other_node) in [(e.u, e.v), (e.v,e.u)]{
                match all_node_info[node as usize].state {
                    0 => all_node_info[node as usize].neigh1=other_node,
                    1 => {
                        all_node_info[node as usize].neigh2=other_node;
                        //Our "node" is now inside the path: we delete all its adjacent edges
                        for pair_id in &mut all_node_info[node as usize].adjacent_pairs{
                            sampled_pairs.delete_pair(*pair_id);
                        }
                        neighborhoods.delete_neighborhood(node);
                        all_node_info[node as usize].adjacent_pairs=Vec::new();   //we clear the vector of adjacent pairs
                        
                    },   
                    2 => panic!("Error. Vertex {node} already in the path"),
                    _ => panic!("Error. How did we end up here??"),
                }
                all_node_info[node as usize].state+=1;
            }

            //Update the endpoints
            let endpu= all_node_info[u as usize].endpoint;
            let endpv= all_node_info[v as usize].endpoint;
            all_node_info[endpu as usize].endpoint=endpv;
            all_node_info[endpv as usize].endpoint=endpu;
        }
    }
}




//calls compute_path_iteration then delete the nodes that are already in some path until there are only two vertices left 
fn compute_path(g:&Graph, all_node_info: &mut Vec<NodeInfo>, nodelist: &mut Vec<Node>, alpha:f32, edgelist: &Vec<Edge>)->f32{
    let mut neighborhoods=MyNeighborhoods::create_neighborhoods(g.n() as usize);
    let mut timer:f32; 
    let mut timetaken:f32=0.0;
    let mut iters:u32=0;

    while nodelist.len()>=3{
        iters+=1;
        println!("Iteration {iters}");
        let now = Instant::now(); // Start time measurement
        compute_path_iteration(&nodelist, g, all_node_info, &mut neighborhoods, alpha, iters, edgelist);
        timer=now.elapsed().as_secs_f32();
        println!("Finished in {timer} sec.");
        timetaken+=timer;
    

        //Remove nodes that are internal in some path from nodelist
        let mut i=0; 
        while i< nodelist.len(){
            if all_node_info[nodelist[i] as usize].state == 2 {
                nodelist.swap_remove(i);
            }
            else{i+=1;}
        }
    }
    timetaken
}








fn main(){
    let args: Vec<String> = env::args().collect();

    //get the graph
    let file_path = &args[1];
    let file_string = fs::read_to_string(file_path)
        .expect("Should have been able to read the file");
    let mut edges: Vec<Edge>= Vec::new();
    for line in file_string.lines() {
        let mut parts = line.split_whitespace(); 
        let u: Node = parts.next().expect("Invalid edge format").parse().expect("Invalid number");
        let v: Node = parts.next().expect("Invalid edge format").parse().expect("Invalid number");
        let edg=Edge{u,v};
        edges.push(edg);
    }
    let g=Graph::from_edges(&mut edges, true);
    println!("n={}, m={}", g.n(), g.m(true));  


    let output_file_name=&args[2];
    let alpha = f32::from_str(&args[3]).unwrap();
    

    /* 
    let first_pair_sampling_path=&args[4];
    let fps_string = fs::read_to_string(first_pair_sampling_path)
    .expect("Should have been able to read the file");
    let mut edgelist: Vec<Edge>= Vec::new();
    for line in fps_string.lines() {
        let mut parts = line.split_whitespace(); 
        let u: Node = parts.next().expect("Invalid edge format").parse().expect("Invalid number");
        let v: Node = parts.next().expect("Invalid edge format").parse().expect("Invalid number");
        let edg=Edge{u,v};
        edgelist.push(edg);
    }*/
    let edgelist: Vec<Edge>= Vec::new();


    //create all_node_info: a vector containing all the nodes of the graph, and some useful info for the algorithm (cf definition of the NodeInfo structure)
    let mut all_node_info: Vec<NodeInfo>=Vec::new();
    for i in 0.. g.n(){
        all_node_info.push(NodeInfo::create_node_info(i));
    }

    //create nodelist: list of nodes we still have to process in compute_path_iteration
    let mut nodelist: Vec<Node>=Vec::new();
    for i in 0.. g.n(){
        nodelist.push(i);
    }


    //Call to the main function
    let timetaken=compute_path(&g, &mut all_node_info, &mut nodelist, alpha, &edgelist);



    //get the results
    assert!(nodelist.len()==2);
    let mut spanning_path: Vec<Node>=vec![0; g.n() as usize];
    let mut all_crossing_numbers=vec![0;g.n() as usize];
    get_results(&all_node_info,&mut spanning_path, &mut all_crossing_numbers, nodelist[0],nodelist[1],&g);
    let max_crossing_number= all_crossing_numbers.iter().max().unwrap_or(&0);

    //In this file, comparison between the degree of each node and its crossing number
    let mut output_file = BufWriter::new(File::create(output_file_name.to_owned()+"_deg_cross.txt").unwrap());
    for i in 0..g.n() {
        writeln!(output_file, "{} {}", g.deg(i), all_crossing_numbers[i as usize]).unwrap();
    }
    output_file.flush().unwrap();


    //in this file the path we obtained
    let mut output_file_path = BufWriter::new(File::create(output_file_name.to_owned()+"_path.txt").unwrap());
    for new_rank in spanning_path {
        writeln!(output_file_path, "{new_rank}").unwrap();
    }
    output_file.flush().unwrap();



    // Some last info
    println!("\t Finished within {:.2?}sec, crossing number: {}", timetaken, max_crossing_number);

    let output_data = OpenOptions::new()
        .append(true)
        .create(true)
        .open("outputdata.txt")
        .unwrap();
    let mut output_data = BufWriter::new(output_data);

    write!(output_data, "{} input size: {} alpha: {} time: {:.2?}sec, crossing nr: {}\n", *file_path, g.n(), alpha, timetaken, max_crossing_number).unwrap();
    output_data.flush().unwrap();
    
}





