#!/usr/bin/env ruby

require 'gnuplot'
require 'set'
require 'byebug'
require 'vose'
require 'pqueue'

class Graph

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
    @degrees = Array.new num_nodes, 0
    parsed_lines.each do |tuple|
      a = tuple[0]
      b = tuple[1]
      @degrees[a] += 1
    end
    @degrees
  end

  def get_num_of_isolated_nodes
    @degrees.count(0)
  end

  def get_max_degree
    @degrees.max
  end

  def get_min_degree
    @degrees.min
  end

  def get_average_degree
    @degrees.inject(:+).fdiv @degrees.length
  end

  def get_density
    num_nodes = @degrees.length
    (2*@degrees.inject(:+))/(num_nodes(num_nodes-1))
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

  def prepare_graph graph_array
    doubled_array = graph_array.dup
    graph_array.each do |line|
      doubled_array << line.reverse
    end
    graph_array = doubled_array.uniq
    graph_array.sort do |x, y|
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
  end

  def prepare_raw_data graph_lines
    parsed_lines = []
    graph_lines.each do |line|
      a, b = line.split(" ").map { |string| string.to_i }
      parsed_lines << [a, b]
    end
    # Delete loops
    parsed_lines.select! do |item|
      item[0] != item[1]
    end
    parsed_lines
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
  end

  def build_graph parsed_lines
    @degrees = get_degrees parsed_lines
    @separators = calculate_separators_in_graph
    @graph = Array.new parsed_lines.length
    parsed_lines.length.times do |index|
      tuple = parsed_lines[index]
      b = tuple[1]
      @graph[index] = b
    end
    @graph
  end

  def calculate_separators_in_graph
    @separators = Array.new @degrees.length
    @separators[0] = 0
    (1..@degrees.length-1).each do |i|
      @separators[i] = @separators[i-1] + @degrees[i-1]
    end
    @separators
  end

  def print_graph_with_separators
    separators = calculate_separators_in_graph.dup
    separators.shift
    separators.map! { |item| item-1 }
    @graph.length.times do |i|
      print @graph[i]
      if separators.include? i
        print "|"
      else
        print ","
      end
    end
    puts ""
  end

  def self.plot_stuff stuff_to_plot, scale='lin', title='', x_label=nil, y_label=nil, style='points'
    Gnuplot.open do |gp|
      Gnuplot::Plot.new(gp) do |plot|

        plot.title title
        plot.xlabel x_label if x_label
        plot.ylabel y_label if y_label

        if stuff_to_plot.is_a? Hash
          x = stuff_to_plot.keys
          y = stuff_to_plot.values
        else
          x = stuff_to_plot[0]
          y = stuff_to_plot[1]
        end
        # if style == 'lines'
        #   second_value = x[1]
        #   x.map! { |thing| thing/second_value }
        # end

        if scale == 'log'
          plot.arbitrary_lines << "set logscale xy"
        end
        plot.arbitrary_lines << "set terminal png size 600,600 enhanced font \"Helvetica,15\""
        plot.arbitrary_lines << "set output \""+title+".png\""

        plot.data << Gnuplot::DataSet.new( [x, y] ) do |ds|
          ds.with = "points" if style == 'points'
          ds.with = "lines" if style == 'lines'
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

  def get_neighbors node_num
    first_separator = @separators[node_num]
    return @graph[first_separator...@graph.length] if node_num+1 == @separators.length
    second_separator = @separators[node_num+1]
    @graph[first_separator...second_separator]
  end

  # Naive approach
  def calculate_clustering_coefficient node_num
    all_neighbors = get_neighbors(node_num).to_set
    # Take care, this number is already two times and should always be even
    all_connections_among_neighbors = 0
    all_neighbors.each do |current_node|
      node_neighbors = get_neighbors(current_node)
      all_connections_among_neighbors += all_neighbors.intersection(node_neighbors).length
    end
    # puts "All connections among neighbors: "+all_connections_among_neighbors.to_s
    # puts "Node degree: "+@degrees[node_num].to_s
    return all_connections_among_neighbors.fdiv(@degrees[node_num]*(@degrees[node_num]-1))
  end

  # Naive approach for clustering coefficient
  def calculate_global_clustering_coefficient
    good_indices = []
    @degrees.length.times do |i|
      if @degrees[i] >= 2
        good_indices << i
      end
    end
    clustering_coefficient_array = (good_indices).to_a.map do |node|
      calculate_clustering_coefficient node
    end
    # puts "All cc #{clustering_coefficient_array.inject(:+)}, degrees.length: #{@degrees.length}"
    clustering_coefficient = clustering_coefficient_array.inject(:+) / @degrees.length
  end

  def check_neighbors_and_degrees
    puts @separators.length
    puts @degrees.length
    @degrees.length.times do |node_num|
      neighbors = get_neighbors node_num
      puts "node_num: "+node_num.to_s+" "+neighbors.length.to_s+" "+@degrees[node_num].to_s
      if neighbors.length != @degrees[node_num]
        return false
      end
    end
    true
  end

  def calculate_triangles_for_each_node
    n = @separators.length
    tr = Array.new n, 0
    n.times do |u|
      neighbors = get_neighbors u, @graph, @separators
      neighbors.each do |v|
        if u < v
          list_U = get_neighbors u, @graph, @separators
          list_V = get_neighbors v, @graph, @separators
          i_u, i_v = 0, 0
          while i_u < @degrees[u] && i_v < @degrees[v] && list_U[i_u] < u && list_V[i_v] < u
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

  def calculate_number_of_triplets_for_each_node
    n = @separators.length
    table = Array.new n, 0
    n.times do |u|
      neighbors = get_neighbors u, @graph, @separators
      neighbors.each do |neighbor|
        table[u] += get_neighbors(neighbor, @graph, @separators).length-1
      end
    end
    table
  end

  def calculate_transitive_ratio
    triangles = calculate_triangles_for_each_node(graph, @separators, @degrees)
    triplets = calculate_number_of_triplets_for_each_node(graph, @separators)
    # puts "Triangles: "+triangles.inject(:+).to_s
    # p triangles
    # puts "Triplets: "+triplets.inject(:+).to_s
    # p triplets
    triangles.inject(:+).fdiv(triplets.inject(:+).fdiv(2))
  end

  def breadth_first_search node
    bfs_output = []
    queue = []
    marked_nodes ||= Set.new
    queue << node
    marked_nodes << node
    while !queue.empty?
      u = queue.shift
      bfs_output << u
      neighbors = get_neighbors u
      neighbors.each do |neighbor|
        if !marked_nodes.include? neighbor
          queue << neighbor
          marked_nodes << neighbor
        end
      end
    end
    bfs_output
  end

  def calculate_betweenness_centrality
    n = @separators.length
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
        get_neighbors(v, @graph, @separators).each do |neighbor|
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

  def get_all_connected_components
    n = @separators.length
    components = []
    already_used_nodes = Set.new
    n.times do |node|
      if already_used_nodes.include? node
        next
      end
      bfs = breadth_first_search node
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

  def generate_hungarian_graph n, m
    graph_array = []
    until graph_array.length == 2*m
      element = [rand(n), rand(n)]
      if element[0] != element[1]
        graph_array << element
        graph_array << [element[1], element[0]]
      end
    end
    graph_array
  end

  def store_graph_in_file file_name, graph_array
    File.open(file_name, 'w') do |f|
      graph_array.each do |edge|
        f.puts edge[0].to_s+" "+edge[1].to_s
      end
    end
  end

  def breadth_first_search_distances node
    # bfs_output = []
    n = @separators.length
    distances = Array.new n, 0
    distances[node] = 0
    queue = []
    marked_nodes ||= Set.new
    queue << node
    marked_nodes << node
    while !queue.empty?
      u = queue.shift
      # bfs_output << u
      neighbors = get_neighbors u
      neighbors.each do |neighbor|
        if !marked_nodes.include? neighbor
          queue << neighbor
          distances[neighbor] = distances[u] + 1
          marked_nodes << neighbor
        end
      end
    end
    distances
  end

  def get_average_distance component
    n = component.length
    global_counter = 0
    n.times do |node|
      distances = breadth_first_search_distances node
      global_counter += distances.inject :+
    end
    global_counter.fdiv(n*(n-1))
  end

  def get_node_from_edge_pos edge_pos
    last_index = 0
    @separators.length.times do |index|
      if @separators[index] > edge_pos
        break
      end
      last_index = index
    end
    # puts "Edge pos #{edge_pos}, start node recovered: #{last_index}"
    last_index
  end

  def generate_fixed_degree_graph degree_list
    @degrees = degree_list
    calculate_separators_in_graph
    num_of_edges = degree_list.inject :+
    puts "number of edges for fixed degree graph: "+num_of_edges.to_s
    if num_of_edges.odd?
      raise "Total number of edges is odd!"
    end
    # puts @separators.inspect
    # puts num_of_edges.inspect
    already_generated_edges = 0
    already_connected = Array.new num_of_edges, false
    who_is_connected_to_whom = (0...degree_list.length).collect { Set.new }
    graph_array = []
    until already_generated_edges >= num_of_edges
      random_start_edge, random_end_edge = nil, nil
      loop do
        random_start_edge = rand(num_of_edges)
        break if !already_connected[random_start_edge]
      end
      loop do
        random_end_edge = rand(num_of_edges)
        break if !already_connected[random_end_edge]
      end
      start_node = get_node_from_edge_pos random_start_edge
      end_node = get_node_from_edge_pos random_end_edge
      if (start_node == end_node || who_is_connected_to_whom[start_node].include?(end_node)) && already_generated_edges+2 == num_of_edges
        graph_array = []
        already_generated_edges = 0
        already_connected = Array.new num_of_edges, false
        who_is_connected_to_whom = (0...degree_list.length).collect { Set.new }
        puts "Already generated #{already_generated_edges}. Restarting graph calculation from the beginning..."
        next
      elsif (start_node == end_node || who_is_connected_to_whom[start_node].include?(end_node)) && already_generated_edges+2 != num_of_edges
        # puts "Already generated #{already_generated_edges}. Chose already connected node..."
        next
      end
      who_is_connected_to_whom[start_node] << end_node
      who_is_connected_to_whom[end_node] << start_node
      already_connected[random_start_edge] = true
      already_connected[random_end_edge] = true
      already_generated_edges += 2
      graph_array << [start_node, end_node]
      # puts graph_array.inspect
      # puts already_connected.inspect
    end
    graph_array
  end

  def get_edge_index_from_nodes start_node, end_node
    # puts 'start_node, end_node '+start_node.inspect+' '+end_node.inspect
    # puts 'first_index, second_index '+@separators[start_node].inspect+' '+get_neighbors(start_node).find_index(end_node).inspect
    second_index = get_neighbors(start_node).find_index(end_node)
    return nil unless second_index
    result = @separators[start_node]+second_index
    # puts "Start node #{start_node}, end node #{end_node}, edge index #{result}"
    result
  end

  def generate_switched_graph num_of_switches, interval=10000
    num_of_edges = @graph.length
    puts "num_of_edges: #{num_of_edges}"
    num_of_switched_ones = 0
    file_name = "switch_clustering_coefficient.txt"
    File.delete file_name if File.exists? file_name
    f = File.open file_name, 'w'
    intermediate_ccs = {}
    until num_of_switches == num_of_switched_ones
      previous_stop = nil
      if num_of_switched_ones % interval == 0 && previous_stop != num_of_switched_ones
        previous_stop = num_of_switched_ones
        avg_clustering_coefficient = calculate_global_clustering_coefficient
        puts "Average clustering coefficient after #{num_of_switched_ones} switches:"
        puts avg_clustering_coefficient.to_s
        intermediate_ccs[num_of_switched_ones] = avg_clustering_coefficient
        f.puts avg_clustering_coefficient.to_s
      end
      first_start_edge = rand(num_of_edges)
      second_start_edge = rand(num_of_edges)
      # puts 'first_start_edge '+first_start_edge.to_s
      # puts 'second_start_edge '+second_start_edge.to_s

      first_start_node = get_node_from_edge_pos first_start_edge
      first_end_node = @graph[first_start_edge]
      second_start_node = get_node_from_edge_pos second_start_edge
      second_end_node = @graph[second_start_edge]

      first_start_node_neighbors = get_neighbors first_start_node
      second_start_node_neighbors = get_neighbors second_start_node
      if first_start_node == second_end_node || second_start_node == first_end_node ||
        first_start_node_neighbors.include?(second_end_node) || second_start_node_neighbors.include?(first_end_node)
        next
      end
      index_first_from_start = first_start_edge
      index_first_from_end = get_edge_index_from_nodes first_end_node, first_start_node
      index_second_from_start = second_start_edge
      index_second_from_end = get_edge_index_from_nodes second_end_node, second_start_node
      # puts "start\nold edge 1: #{first_start_node}, #{first_end_node}"
      # puts "old edge 2: #{second_start_node}, #{second_end_node}"
      # puts "new edge 1: #{first_start_node}, #{second_end_node}"
      # puts "new edge 2: #{second_start_node}, #{first_end_node}"
      # puts "start\nold value #{@graph[index_first_from_start]}, new value #{second_end_node}"
      @graph[index_first_from_start] = second_end_node
      # puts "old value #{@graph[index_first_from_end]}, new value #{second_start_node}"
      @graph[index_first_from_end] = second_start_node
      # puts "old value #{@graph[index_second_from_end]}, new value #{first_start_node}"
      @graph[index_second_from_end] = first_start_node
      # puts "old value #{@graph[index_second_from_start]}, new value #{first_end_node}"
      @graph[index_second_from_start] = first_end_node
      num_of_switched_ones += 1
    end
    return intermediate_ccs
  end

  def insert_edge start_node, end_node, revert=true
    # puts "insert edge revert=#{revert}"
    neighbors = get_neighbors start_node
    last_neighbor_index = 0
    neighbors.length.times do |current_node_index|
      if neighbors[current_node_index] > end_node
        break
      end
      last_neighbor_index = current_node_index
    end
    edge_index = @separators[start_node]+last_neighbor_index
    @degrees[start_node] += 1
    @graph.insert edge_index, end_node
    @separators.length.times do |separator_index|
      if separator_index > start_node
        @separators[separator_index] += 1
      end
    end
    if revert
      insert_edge end_node, start_node, false
    end
    # byebug
  end

  def insert_node
    # puts "insert node"
    @separators += [@graph.length]
    @degrees += [0]
  end

  def remove_edge start_node, end_node, revert=true
    @degrees[start_node] -= 1
    edge_index = get_edge_index_from_nodes start_node, end_node
    new_separators = []
    @separators.each do |separator|
      if separator > edge_index
        separator -= 1
      end
      new_separators << separator
    end
    @separators = new_separators
    @graph.delete_at edge_index
    if revert
      remove_edge end_node, start_node, false
    end
  end

  def generate_scale_free_graph n, k
    g_new = Graph.new
    g_new.instance_variable_set :@graph, @graph
    g_new.instance_variable_set :@separators, @separators
    g_new.instance_variable_set :@degrees, @degrees
    i = @separators.length
    while i <= n do
      # puts "i: "+i.to_s+", n: "+n.to_s
      sum_of_degrees = g_new.instance_variable_get(:@degrees).inject :+
      # Why should k be in the following formula?
      probabilities = g_new.instance_variable_get(:@degrees).map { |item| item.fdiv sum_of_degrees }
      # puts "degrees: #{g_new.instance_variable_get(:@degrees)}"
      # puts "probabilities: #{probabilities.inspect}, sum: #{probabilities.inject :+}"
      already_chosen_values = []
      vose = Vose::AliasMethod.new probabilities
      until already_chosen_values.length == k
        choice = vose.next
        next if already_chosen_values.include? choice
        already_chosen_values << choice
      end
      g_new.insert_node
      already_chosen_values.each do |end_node|
        g_new.insert_edge i, end_node
      end
      i += 1
    end
    g_new
  end

# Inspired by https://github.com/Florent2/wsmodel/blob/master/lib/wsmodel/network.rb
  def generate_watts_strogatz num, avg_degree
    array = []
    num.times do |node|
      (avg_degree/2).times do |i|
        array << [node, (node+i+1) % num]
      end
    end
    array
  end

  def switch_watts_strogatz avg_degree, beta
    num_nodes = @separators.length
    # puts num_nodes
    (avg_degree/2).times do |i|
      num_nodes.times do |node|
        # puts "i: #{i}, node: #{node}"
        if rand < beta
          neighbor = (node + i + 1) % num_nodes
          # puts "#{node} #{neighbor}"
          remove_edge node, neighbor
          unlinkable_nodes = [node, neighbor] + get_neighbors(node)
          new_neighbor = ((0...num_nodes).to_a - unlinkable_nodes).sample
          insert_edge node, new_neighbor
        end
      end
    end
  end

  def write_to_file file_name
    File.delete(file_name) if File.exists? file_name
    f = File.open file_name, 'w'
    @separators.length.times do |node|
      neighbors = get_neighbors node
      neighbors.each do |neighbor|
        f.puts node.to_s+' '+neighbor.to_s
      end
    end
  end

  def check_graph_sanity
    @separators.length.times do |node|
      neighbors = get_neighbors node
      if neighbors.uniq.length != neighbors.length
        raise "Multiple edge found for node #{node} with neighbors #{neighbors}"
      end
      if neighbors.include? node
        raise "Node #{node} has a loop!"
      end
    end
    puts "Graph is sane!"
  end

  def get_precision_and_recall_data predictions, edges_to_guess
    edges_to_guess = Set.new(edges_to_guess)
    # puts "Rest of graph: #{edges_to_guess}"

    true_positives = Array.new predictions.length, nil
    false_positives = Array.new predictions.length, nil
    false_negatives = Array.new predictions.length, nil
    previous_index = nil
    predictions.length.times do |i|
      tp_zero_or_one = 0
      # byebug
      if edges_to_guess.include?(predictions[i]) ||
        edges_to_guess.include?(predictions[i].reverse)
        # puts "Wohooooo!"
        tp_zero_or_one = 1
      end
      previous_tp_value = 0
      previous_tp_value = true_positives[previous_index] if previous_index
      true_positives[i] = previous_tp_value + tp_zero_or_one

      previous_fp_value = 0
      previous_fp_value = false_positives[previous_index] if previous_index
      false_positives[i] = previous_fp_value + 1-tp_zero_or_one

      false_negatives[i] = edges_to_guess.length - true_positives[i]
      previous_index = i
    end
    return true_positives, false_positives, false_negatives
  end

  def calculate_precision_and_recall tp, fp, fn
    num_of_predictions = tp.length
    pr = []
    rc = []
    num_of_predictions.times do |i|
      pr[i] = tp[i].fdiv(tp[i]+fp[i])
      rc[i] = tp[i].fdiv(tp[i]+fn[i])
    end
    return pr, rc
  end

  def strip_graph_of_all_edges
    sample = Graph.new
    sample.instance_variable_set(:@degrees, @degrees.map{ |degree| 0 })
    sample.instance_variable_set(:@separators, @separators.map{ |separator| 0 })
    sample.instance_variable_set(:@graph, [])
    sample
  end

  def write_to_res_file start_node, end_node, num_of_tests_made
    @res_file.puts "(#{num_of_tests_made},#{start_node},#{end_node})"
    # puts "(#{num_of_tests_made},#{start_node},#{end_node})"
  end

  # To do: Implement efficiency of pure random strategy
  def analyse reference_graph, num_of_tests_made
    num_of_nodes = reference_graph.instance_variable_get(:@degrees).length
    num_of_edges = reference_graph.instance_variable_get(:@graph).length/2
    possible_edges = num_of_nodes*(num_of_nodes-1)/2
    # puts "More tests made than there are possible edges!" if num_of_tests_made > possible_edges
    efficiency_worst = 0
    num_of_tests_made.times do |i|
      thing_to_add = [[i+1-(possible_edges-num_of_edges), num_of_edges].min, 0].max
      # puts thing_to_add
      efficiency_worst += thing_to_add
    end
    efficiency_best = 0
    num_of_tests_made.times do |i|
      thing_to_add = [i+1, num_of_edges].min
      # puts "thing_to_add: #{thing_to_add}"
      # puts "i+1: #{i+1}"
      # puts "num_of_edges: #{num_of_edges}"
      efficiency_best += thing_to_add
    end
    efficiency_random = 0
    number_of_discovered_edges = 0
    num_of_tests_made.times do |i|
      number_of_discovered_edges += (num_of_edges-number_of_discovered_edges).fdiv(possible_edges-i)
      # puts "thing_to_add: #{thing_to_add}"
      # puts "i+1: #{i+1}"
      # puts "num_of_edges: #{num_of_edges}"
      efficiency_random += number_of_discovered_edges
    end
    return efficiency_worst, efficiency_best, efficiency_random
  end

  def calculate_efficiency reference_graph, num_of_tests_made
    num_of_nodes = reference_graph.instance_variable_get(:@degrees).length
    num_of_edges = reference_graph.instance_variable_get(:@graph).length/2
    lines = IO.readlines(@res_file)
    num_of_found_edges = lines.length
    efficiency_worst, efficiency_best = analyse reference_graph, num_of_tests_made
    efficiency = 0
    current_line = 0
    current_thing_to_add = 0
    num_of_tests_made.times do |i|
      i += 1
      matched_stuff = /\d+/.match(lines[current_line])
      next unless matched_stuff
      filtered_number = matched_stuff[0].to_i
      # puts filtered_number
      if i == filtered_number
        current_thing_to_add = filtered_number
        current_line += 1
        efficiency += current_thing_to_add
      end
    end
    return efficiency
  end

  def random_strategy reference_graph, num_of_tests
    num_of_nodes = reference_graph.instance_variable_get(:@degrees).length
    tested_pairs = Set.new
    num_of_tests.times do |test_num|
      start_node, end_node = nil, nil
      loop do
        start_node = rand num_of_nodes
        end_node = rand num_of_nodes
        break unless tested_pairs.include?([start_node, end_node]) || tested_pairs.include?([end_node, start_node])
      end
      edge_exists = reference_graph.get_edge_index_from_nodes start_node, end_node
      # puts "Edge exists: #{edge_exists}"
      tested_pairs << [start_node, end_node]
      if edge_exists
        insert_edge start_node, end_node
        write_to_res_file start_node, end_node, test_num
        # puts "Edges in the graph: #{@graph.length/2}"
      end
    end
    return tested_pairs
  end

  # To do
  # .map.with_index.sort.map(&:last)
  # take tested pairs from random to complete
  def complete_strategy reference_graph, tested_pairs=Set.new
    num_of_nodes = reference_graph.instance_variable_get(:@degrees).length
    tested_pairs = Set.new
    nodes_with_degree = []
    queue = PQueue.new(nodes_with_degree) { |a,b| @degrees[a] > @degrees[b] }
    # queue = @degrees.map.with_index.sort.reverse.map(&:last)[-@degrees.count(0),-1]
    until queue.empty?
      start_node = queue.shift
      num_of_nodes.times do |end_node|
        next if tested_pairs.include?([start_node, end_node]) || tested_pairs.include?([end_node, start_node])
        edge_exists = reference_graph.get_edge_index_from_nodes start_node, end_node
        # puts "Edge exists: #{edge_exists}"
        tested_pairs << [start_node, end_node]
        if edge_exists && @degrees[end_node] == 0
          insert_edge start_node, end_node
          queue << end_node
          queue.send :sort!
          write_to_res_file start_node, end_node, test_num
          # puts "Edges in the graph: #{@graph.length/2}"
        end
      end
    end
    # byebug
  end

  def extract_data_for_plotting
    x_axis = []
    y_axis = []
    found_edges_counter = 0
    # byebug
    @res_file.seek 0
    @res_file.each do |line|
      # puts "spam and eggs"
      found_edges_counter += 1
      test_num = (/\d+/.match line)[0].to_i
      x_axis << test_num
      y_axis << found_edges_counter
    end
    return x_axis, y_axis
  end

end

g_original = Graph.new
prepared_array = g_original.prepare_data File.open("flickr-test", "r").each_line
g_original.build_graph prepared_array
puts "Number of nodes:"
puts g_original.instance_variable_get(:@separators).length.to_s
puts "Number of edges:"
puts g_original.instance_variable_get(:@graph).length.to_s

# sample = g_original.strip_graph_of_all_edges
# num_of_tests_made = 50000
# sample.instance_variable_set(:@res_file, File.open("random_results.txt", 'w+'))
# sample.random_strategy g_original, num_of_tests_made
# puts "Number of found edges:\t#{sample.instance_variable_get(:@graph).length/2}"
# stuff_to_plot = sample.extract_data_for_plotting
# stuff_to_plot[0].map! { |thing| thing.fdiv(1000).round() }
#
# efficiency_worst, efficiency_best, efficiency_random = sample.analyse g_original, num_of_tests_made
# puts "Worst efficiency:\t#{efficiency_worst}"
# puts "Best efficiency:\t#{efficiency_best}"
# puts "Random Efficiency:\t#{efficiency_random}"
# efficiency = sample.calculate_efficiency g_original, num_of_tests_made
# puts "Efficiency:\t\t#{efficiency}"
# normalized_efficiency = (efficiency - efficiency_worst).fdiv(efficiency_best - efficiency_worst)
# normalized_efficiency_random = (efficiency_random - efficiency_worst).fdiv(efficiency_best - efficiency_worst)
# puts "Normalized efficiency:\t#{normalized_efficiency}"
# puts "Relative efficiency:\t#{normalized_efficiency/normalized_efficiency_random}"
# # Graph::plot_stuff stuff_to_plot, scale='lin', title='Random', x_label="number of tests in thousands", y_label="found edges", style='lines'

sample = g_original.strip_graph_of_all_edges
num_of_tests_made = 50000
sample.instance_variable_set(:@res_file, File.open("complete_results.txt", 'w+'))
tested_pairs = sample.random_strategy g_original, num_of_tests_made
sample.complete_strategy g_original, tested_pairs
puts "Number of found edges:\t#{sample.instance_variable_get(:@graph).length/2}"
stuff_to_plot = sample.extract_data_for_plotting
stuff_to_plot[0].map! { |thing| thing.fdiv(1000).round() }

efficiency_worst, efficiency_best, efficiency_random = sample.analyse g_original, num_of_tests_made
puts "Worst efficiency:\t#{efficiency_worst}"
puts "Best efficiency:\t#{efficiency_best}"
puts "Random Efficiency:\t#{efficiency_random}"
efficiency = sample.calculate_efficiency g_original, num_of_tests_made
puts "Efficiency:\t\t#{efficiency}"
normalized_efficiency = (efficiency - efficiency_worst).fdiv(efficiency_best - efficiency_worst)
normalized_efficiency_random = (efficiency_random - efficiency_worst).fdiv(efficiency_best - efficiency_worst)
puts "Normalized efficiency:\t#{normalized_efficiency}"
puts "Relative efficiency:\t#{normalized_efficiency/normalized_efficiency_random}"
# Graph::plot_stuff stuff_to_plot, scale='lin', title='Random', x_label="number of tests in thousands", y_label="found edges", style='lines'
