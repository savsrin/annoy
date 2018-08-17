struct Blosum {
  template<typename S, typename T>
  struct ANNOY_NODE_ATTRIBUTE Node {
    S n_descendants; //type int 
    S children[2]; //array of two ints/longs
    T v[1]; // array of one float
  };

  static const size_t max_iterations = 20;
  
  static const int num_amino_acids = 24; 

  static const float scores62 [num_amino_acids] [num_amino_acids] = {
    {4, -1, -2, -2, 0, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -3, -2, 0, -2, -1, 0, -4},
    {-1, 5, 0, -2, -3, 1, 0, -2, 0, -3, -2, 2, -1, -3, -2, -1, -1, -3, -2, -3, -1, 0, -1, -4},
    {-2, 0, 6, 1, -3, 0, 0, 0, 1, -3, -3, 0, -2, -3, -2, 1, 0, -4, -2, -3, 3, 0, -1, -4}, 
    {-2, -2, 1, 6, -3, 0, 2, -1, -1, -3, -4, -1, -3, -3, -1, 0, -1, -4, -3, -3, 4, 1, -1, -4},  
    {0, -3, -3, -3, 9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -3, -2, -4},
    {-1, 1, 0, 0, -3, 5, 2, -2, 0, -3, -2, 1, 0, -3, -1, 0, -1, -2, -1, -2, 0, 3, -1, -4}, 
    {-1, 0, 0, 2, -4, 2, 5, -2, 0, -3, -3, 1, -2, -3, -1, 0, -1, -3, -2, -2, 1, 4, -1, -4},
    {0, -2, 0, -1, -3, -2, -2, 6, -2, -4, -4, -2, -3, -3, -2, 0, -2, -2, -3, -3, -1, -2, -1, -4},
    {-2, 0, 1, -1, -3, 0, 0, -2, 8, -3, -3, -1, -2, -1, -2, -1, -2, -2, 2, -3, 0, 0, -1, -4},
    {-1, -3, -3, -3, -1, -3, -3, -4, -3, 4, 2, -3, 1, 0, -3, -2, -1, -3, -1, 3, -3, -3, -1, -4},
    {-1, -2, -3, -4, -1, -2, -3, -4, -3,  2, 4, -2, 2, 0, -3, -2, -1, -2, -1, 1, -4, -3, -1, -4},
    {-1, 2, 0, -1, -3, 1, 1, -2, -1, -3, -2, 5, -1, -3, -1, 0, -1, -3, -2, -2, 0, 1, -1 -4},
    {-1, -1, -2, -3, -1, 0, -2, -3, -2, 1, 2, -1, 5, 0, -2, -1, -1, -1, -1, 1, -3, -1, -1, -4},
    {-2, -3, -3, -3, -2, -3, -3, -3, -1, 0, 0, -3, 0, 6, -4, -2, -2, 1, 3, -1, -3, -3, -1, -4},
    {-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4, 7, -1, -1, -4, -3, -2, -2, -1, -2, -4},
    {1, -1, 1, 0, -1, 0, 0, -1, -2, -2, 0, -1, -2, -1, 4, 1, -3, -2, -2, 0, 0, 0, -4},
    {0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1, 5, -2, -2, 0, -1, -1, 0, -4},
    {-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1, -4, -3, -2, 11, 2, -3, -4, -3, -2, -4},
    {-2, -2, -2, -3, -2, -1, -2, -3, 2, -1, -1, -2, -1, 3, -3, -2, -2, 2, 7, -1, -3, -2, -1, -4},
    {0, -3, -3, -3, -1, -2, -2, -3, -3, 3, 1, -2, 1, -1, -2, -2, 0, -3, -1, 4, -3, -2, -1, -4},
    {-2, -1, 3, 4, -3, 0, 1, -1, 0, -3, -4, 0, -3, -3, -2, 0, -1, -4, -3, -3, 4, 1, -1, -4},
    {-1, 0, 0, 1, -3, 3, 4, -2, 0, -3, -3, 1, -1, -3, -1, 0, -1, -3, -2, -2, 1, 4, -1, -4},
    {0, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2, 0, 0, -2, -1, -1, -1, -1, -1, -4},
    {0, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2, 0, 0, -2, -1, -1, -1, -1, -1, -4}
   
  };

  template<typename T>
  static inline T pq_distance(T distance, T margin, int child_nr) {
    return distance - (margin != (unsigned int) child_nr);
  }

  template<typename T>
  static inline T pq_initial_value() {
    return numeric_limits<T>::max();
  }
  template<typename S, typename T>
  static inline T distance(const Node<S, T>* x, const Node<S, T>* y, int f) {
    size_t dist = 0;
    for (int i = 0; i < f; i++) {
      dist += scores62[x->v[i]][y->v[i]];
    }
    return dist;
  }
  template<typename S, typename T>
  static inline bool margin(const Node<S, T>* n, const T* y, int f) {
    static const size_t n_bits = sizeof(T) * 8;
    T chunk = n->v[0] / n_bits;
    return (y[chunk] & (static_cast<T>(1) << (n_bits - 1 - (n->v[0] % n_bits)))) != 0;
  }
  template<typename S, typename T, typename Random>
  static inline bool side(const Node<S, T>* n, const T* y, int f, Random& random) {
    return margin(n, y, f);
  }
  template<typename S, typename T, typename Random>
  static inline void create_split(const vector<Node<S, T>*>& nodes, int f, size_t s, Random& random, Node<S, T>* n) {
    size_t cur_size = 0;
    size_t i = 0;
    int dim = f * 8 * sizeof(T);
    for (; i < max_iterations; i++) {
      // choose random position to split at
      n->v[0] = random.index(dim);
      cur_size = 0;
      for (typename vector<Node<S, T>*>::const_iterator it = nodes.begin(); it != nodes.end(); ++it) {
        if (margin(n, (*it)->v, f)) {
          cur_size++;
        }
      }
      if (cur_size > 0 && cur_size < nodes.size()) {
        break;
      }
    }
    // brute-force search for splitting coordinate
    if (i == max_iterations) {
      int j = 0;
      for (; j < dim; j++) {
        n->v[0] = j;
        cur_size = 0;
  for (typename vector<Node<S, T>*>::const_iterator it = nodes.begin(); it != nodes.end(); ++it) {
          if (margin(n, (*it)->v, f)) {
            cur_size++;
          }
        }
        if (cur_size > 0 && cur_size < nodes.size()) {
          break;
        }
      }
    }
  }
  template<typename T>
  static inline T normalized_distance(T distance) {
    return distance;
  }
  template<typename S, typename T>
  static inline void init_node(Node<S, T>* n, int f) {
  }
  static const char* name() {
    return "blosum";
  }
};

