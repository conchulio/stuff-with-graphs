#!/usr/bin/env ruby

require 'gnuplot'
require 'set'

def get_num_of_nodes parsed_lines
  unique_things = []
  parsed_lines.each do |tuple|
    a = tuple[0]
    b = tuple[1]
    unique_things << a
    unique_things << b
  end
  unique_things.max + 1
end

def get_degrees parsed_lines
  num_nodes = get_num_of_nodes parsed_lines
  # degrees = [0]*num_nodes
  degrees = Array.new num_nodes, 0
  parsed_lines.each do |tuple|
    a = tuple[0]
    b = tuple[1]
    degrees[a] += 1
  end
  degrees
end

def get_num_of_isolated_nodes degrees
  degrees.count(0)
end

def get_max_degree degrees
  degrees.max
end

def get_min_degree degrees
  degrees.min
end

def get_average_degree degrees
  degrees.inject(:+)/degrees.length
end

def get_density degrees
  num_nodes = degrees.length
  (2*degrees.inject(:+))/(num_nodes(num_nodes-1))
end

def get_distribution things
  distribution = {}
  things.each do |thing|
    if !distribution.has_key? thing
      distribution[thing] = 0
    end
    distribution[thing] += 1
  end
  distribution
end

def get_degree_distribution degrees
  get_distribution degrees
end

def prepare_data graph_lines
  parsed_lines = []
  graph_lines.each do |line|
    a, b = line.split(" ").map { |string| string.to_i }
    parsed_lines << [a, b]
  end
  # Delete loops
  parsed_lines.select! do |item|
    item[0] != item[1]
  end
  new_parsed_lines = []
  parsed_lines.each do |item|
    new_parsed_lines << item
    new_parsed_lines << [item[1], item[0]]
  end
  parsed_lines = new_parsed_lines.uniq
  parsed_lines.sort! do |x, y|
    if x[0] < y[0]
      -1
    elsif x[0] == y[0]
      if x[1] < y[1]
        -1
      elsif x[1] == y[1]
        0
      else
        +1
      end
    else
      +1
    end
  end
  doubled_lines = parsed_lines.dup
  parsed_lines.each do |line|
    doubled_lines << line.reverse
  end
  doubled_lines.uniq!
  doubled_lines
end

def build_graph parsed_lines
  degrees = get_degrees parsed_lines
  separators = calculate_separators_in_graph degrees
  graph = Array.new parsed_lines.length
  (0..parsed_lines.length-1).each do |index|
    tuple = parsed_lines[index]
    b = tuple[1]
    graph[index] = b
  end
  graph
end

def calculate_separators_in_graph degrees
  separators = Array.new degrees.length
  separators[0] = 0
  (1..degrees.length-1).each do |i|
    separators[i] = separators[i-1] + degrees[i-1]
  end
  separators
end

def print_graph_with_separators degrees, graph
  separators = calculate_separators_in_graph degrees
  separators.shift
  separators.map! { |item| item-1 }
  (0..graph.length-1).each do |i|
    print graph[i]
    if separators.include? i
      print "|"
    else
      print ","
    end
  end
  puts ""
end

def plot_stuff degrees_dist, scale='lin', title=''
  Gnuplot.open do |gp|
    Gnuplot::Plot.new(gp) do |plot|

      plot.title title
      # plot.xlabel "x"
      # plot.ylabel "y"

      x = degrees_dist.keys
      y = degrees_dist.values

      if scale == 'log'
        plot.arbitrary_lines << "set logscale xy"
      end
      plot.arbitrary_lines << "set terminal png size 600,600 enhanced font \"Helvetica,15\""
      plot.arbitrary_lines << "set output \""+title+".png\""

      plot.data << Gnuplot::DataSet.new( [x, y] ) do |ds|
        ds.with = "points"
        ds.notitle
      end
    end
  end
end

def get_cumulative_distribution degrees_dist
  cumulative_distribution = {}
  value_from_previous = 0
  degrees_dist.keys.sort.each do |degree|
    cumulative_distribution[degree] = value_from_previous + degrees_dist[degree]
    value_from_previous = cumulative_distribution[degree]
  end
  cumulative_distribution
end

def get_neighbors node_num, graph, separators
  if node_num+1 == separators.length
    # puts "I'm in this case, separator: "+separators[node_num].to_s+", second separator: "+separators.length.to_s
    all_neighbors = graph[separators[node_num]..graph.length-1]
  elsif separators[node_num] == separators[node_num+1]
    all_neighbors = []
  else
    all_neighbors = graph[separators[node_num]..separators[node_num+1]-1]
  end
end

# Naive approach
def calculate_clustering_coefficient node_num, degrees, graph, separators
  all_neighbors = get_neighbors(node_num, graph, separators).to_set
  # Take care, this number is already two times and should always be even
  all_connections_among_neighbors = 0
  all_neighbors.each do |current_node|
    node_neighbors = get_neighbors(current_node, graph, separators)
    all_connections_among_neighbors += all_neighbors.intersection(node_neighbors).length
  end
  # puts "All connections among neighbors: "+all_connections_among_neighbors.to_s
  # puts "Node degree: "+degrees[node_num].to_s
  return all_connections_among_neighbors.fdiv(degrees[node_num]*(degrees[node_num]-1))
end

# Naive approach for clustering coefficient
def calculate_global_clustering_coefficient degrees, graph, separators
  good_indices = []
  (0..degrees.length-1).each do |i|
    if degrees[i] >= 2
      good_indices << i
    end
  end
  clustering_coefficient_array = (good_indices).to_a.map do |node|
    calculate_clustering_coefficient(node, degrees, graph, separators)
  end
  clustering_coefficient = clustering_coefficient_array.instance_eval { reduce(:+) / degrees.length }
end

def check_neighbors_and_degrees graph, separators, degrees
  puts separators.length
  puts degrees.length
  (0..degrees.length-1).each do |node_num|
    neighbors = get_neighbors(node_num, graph, separators)
    puts "node_num: "+node_num.to_s+" "+neighbors.length.to_s+" "+degrees[node_num].to_s
    if neighbors.length != degrees[node_num]
      return false
    end
  end
  true
end

def calculate_triangles_for_each_node graph, separators, degrees
  n = separators.length
  tr = Array.new n, 0
  (0..n-1).each do |u|
    neighbors = get_neighbors u, graph, separators
    neighbors.each do |v|
      if u < v
        list_U = get_neighbors u, graph, separators
        list_V = get_neighbors v, graph, separators
        i_u, i_v = 0, 0
        while i_u < degrees[u] && i_v < degrees[v] && list_U[i_u] < u && list_V[i_v] < u
          if list_U[i_u] < list_V[i_v]
            i_u += 1
          else
            if list_U[i_u] > list_V[i_v]
              i_v += 1
            else
              tr[u] += 1
              tr[v] += 1
              tr[list_U[i_u]] += 1
              i_u += 1
              i_v += 1
            end
          end
        end
      end
    end
  end
  tr
end

def calculate_number_of_triplets_for_each_node graph, separators
  n = separators.length
  table = Array.new n, 0
  (0..n-1).each do |u|
    neighbors = get_neighbors u, graph, separators
    neighbors.each do |neighbor|
      table[u] += get_neighbors(neighbor, graph, separators).length-1
    end
  end
  table
end

def calculate_transitive_ratio graph, separators, degrees
  triangles = calculate_triangles_for_each_node(graph, separators, degrees)
  triplets = calculate_number_of_triplets_for_each_node(graph, separators)
  # puts "Triangles: "+triangles.inject(:+).to_s
  # p triangles
  # puts "Triplets: "+triplets.inject(:+).to_s
  # p triplets
  triangles.inject(:+).fdiv(triplets.inject(:+).fdiv(2))
end

def breadth_first_search node, graph, separators
  bfs_output = []
  queue = []
  marked_nodes ||= Set.new
  queue << node
  marked_nodes << node
  while !queue.empty?
    u = queue.shift
    bfs_output << u
    neighbors = get_neighbors u, graph, separators
    neighbors.each do |neighbor|
      if !marked_nodes.include? neighbor
        queue << neighbor
        marked_nodes << neighbor
      end
    end
  end
  bfs_output
end

# def calculate_betweenness_centrality graph, separators
#   n = separators.length
#   cb = Array.new n, 0
#   n.times do |s|
#   	cs = []
#   	cp = (0...n).collect { [] }
#   	sigma = Array.new n, 0
#   	sigma[s] = 1
#   	d = Array.new n, -1
#   	d[s] = 0
#   	cq = []
#   	cq << s # enqueue
#   	while !cq.empty?
#   		v = cq.shift
#   		cs.unshift v
#   		get_neighbors(v, graph, separators).each do |w|
#   			if d[w] < 0
#   				cq << w
#   				d[w] = d[v] + 1
#   			end
#   			if d[w] == d[v] + 1
#   				sigma[w] += sigma[v]
#   				cp[w] << v
#   			end
#   		end
#   	end
#   	delta = Array.new n, 0
#   	cs.each do |w|
#   		cp[w].each do |v|
#   			delta[v] += (sigma[v].to_f / sigma[w]) * (1+delta[w])
#   			cb[w] += delta[w] if w != s
#   		end
#   	end
#   end
#   cb
# end

def calculate_betweenness_centrality graph, separators
  n = separators.length
  c_b = Array.new n, 0
  n.times do |s|
    stack = []
    list_p = (0...n).collect { [] }
    sigma = Array.new n, 0
    sigma[s] = 1
    d = Array.new n, -1
    d[s] = 0
    q = []
    q << s
    while !q.empty?
      v = q.shift
      stack << v
      get_neighbors(v, graph, separators).each do |neighbor|
        if d[neighbor] < 0
          q << neighbor
          d[neighbor] = d[v] + 1
        end
        if d[neighbor] == d[v] + 1
          sigma[neighbor] = sigma[neighbor] + sigma[v]
          list_p[neighbor] << v
        end
      end
    end
    delta = Array.new n, 0
    while !stack.empty?
      w = stack.pop
      list_p[w].each do |some_node|
        delta[some_node] = delta[some_node] + (sigma[some_node].fdiv(sigma[w]))*(1+delta[w])
      end
      if w != s
        c_b[w] = c_b[w] + delta[w]
      end
    end
  end
  c_b
end

def all_connected_components graph, separators
  n = separators.length
  components = []
  already_used_nodes = Set.new
  (0..n-1).each do |node|
    if already_used_nodes.include? node
      next
    end
    bfs = breadth_first_search node, graph, separators
    already_used_nodes += bfs
    components << bfs
  end
  components.sort do |x, y|
    if x.length < y.length
      1
    elsif x.length == y.length
      0
    else
      -1
    end
  end
end

parsed_lines = prepare_data ARGF
num_nodes = get_num_of_nodes parsed_lines

puts "Exercise 2, number of nodes:"
puts num_nodes.inspect
File.open('graphe.n','w') do |s|
  s.puts num_nodes.inspect
end
degrees = get_degrees parsed_lines

puts "Exercise 3, node degrees:"
puts degrees.inspect
File.open('graphe.deg','w') do |s|
  s.puts degrees.inspect
end

puts "Exercise 5, number of isolated nodes:"
puts get_num_of_isolated_nodes(degrees).to_s

degrees_dist = get_degree_distribution degrees
puts "Exercise 6, degree distribution:"
printable_degree_dist = degrees_dist.to_a.sort do |x, y|
  if x[0] < y[0]
    -1
  elsif x[0] == y[0]
    if x[1] < y[1]
      -1
    elsif x[1] == y[1]
      0
    else
      +1
    end
  else
    +1
  end
end
printable_degree_dist.each do |tuple|
  puts tuple[0].to_s+" "+tuple[1].to_s
end
File.open('graphe.dist','w') do |s|
  printable_degree_dist.each do |tuple|
    s.puts tuple[0].to_s+" "+tuple[1].to_s
  end
end

# plot_stuff degrees_dist, 'log', 'degree distribution'
# cumulative_distribution = get_cumulative_distribution degrees_dist
# plot_stuff cumulative_distribution, 'log', 'cumulative distribution'

graph = build_graph parsed_lines
separators = calculate_separators_in_graph degrees

global_clustering_coefficient = calculate_global_clustering_coefficient degrees, graph, separators
puts "Exercise 10, global clustering coefficient:"
puts global_clustering_coefficient.to_s
transitive_ratio = calculate_transitive_ratio graph, separators, degrees
puts "Exercise 10, transitive_ratio:"
puts transitive_ratio.to_s

puts "Exercise 12, distribution of connected components:"
components = all_connected_components(graph, separators)
components_distribution = get_distribution(components.map { |item| item.length })
printable_components = components_distribution.to_a.sort do |x, y|
  if x[0] < y[0]
    -1
  elsif x[0] == y[0]
    if x[1] < y[1]
      -1
    elsif x[1] == y[1]
      0
    else
      +1
    end
  else
    +1
  end
end
# plot_stuff components_distribution, 'log', 'components distribution'
printable_components.each do |tuple|
  puts tuple[0].to_s+" "+tuple[1].to_s
end
puts "Exercise 13, principal connected component:"
puts components[0].inspect

# bfs_output = breadth_first_search node_for_bfs, graph, separators
# puts "BFS from node: "+node_for_bfs.to_s+", "+bfs_output.inspect

# components = all_connected_components graph, separators
# puts "Size of components: "+components.map{ |c| c.length }.inspect
# puts "Sum of all components: "+components.map{ |c| c.length }.inject(:+).to_s

# puts "Neighbor retrieval works?: "+check_neighbors_and_degrees(graph, separators, degrees).to_s

# (0..num_nodes-1).each do |node|
#   puts "Node: "+node.to_s+", neighbors: "+get_neighbors(node, graph, separators).inspect
# end

betweenness_centrality = calculate_betweenness_centrality(graph, separators)
puts "Exercise 16, beetweenness centrality:"
puts betweenness_centrality.inspect
bc_dist = get_distribution(betweenness_centrality)
printable_bc_dist = bc_dist.to_a.sort do |x, y|
  if x[0] < y[0]
    -1
  elsif x[0] == y[0]
    if x[1] < y[1]
      -1
    elsif x[1] == y[1]
      0
    else
      +1
    end
  else
    +1
  end
end
# plot_stuff components_distribution, 'log', 'components distribution'
puts "Exercise 16, distribution of betweenness centralities:"
printable_bc_dist.each do |tuple|
  puts tuple[0].to_s+" "+tuple[1].to_s
end
cumulative_bc_dist = get_cumulative_distribution bc_dist
puts "Exercise 16, cumulative distribution of betweenness centralities:"
printable_cumulative_bc_dist = cumulative_bc_dist.to_a.sort do |x, y|
  if x[0] < y[0]
    -1
  elsif x[0] == y[0]
    if x[1] < y[1]
      -1
    elsif x[1] == y[1]
      0
    else
      +1
    end
  else
    +1
  end
end
printable_cumulative_bc_dist.each do |tuple|
  puts tuple[0].to_s+" "+tuple[1].to_s
end
plot_stuff bc_dist, 'log', 'betweenness centrality distribution'
plot_stuff cumulative_bc_dist, 'log', 'betweenness centrality cumulative distribution'
