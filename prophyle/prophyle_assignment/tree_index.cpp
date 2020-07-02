#include "tree_index.h"
#include "word_splitter.h"
#include "knhx.h"
#include <sstream>
#include <fstream>
#include <cstring>

TreeIndex::TreeIndex(const std::string& tree_filename) {
  std::ifstream t(tree_filename);
  std::stringstream buffer;
  buffer << t.rdbuf();
  char s[buffer.str().length()+1];
  strcpy(s, buffer.str().c_str());
  int error;
  first_node_ = kn_parse(s, &nodes_count_, &error);
  bool root_found = false;
  for (int32_t i = 0; i < nodes_count_; ++i) {
    knhx1_t* possible_root = first_node_ + i;
    if (possible_root->parent == -1) {
      root_ = possible_root;
      root_found = true;
      break;
    }
  }
  if (!root_found) {
    std::cerr << "root not found" << std::endl;
    exit(1);
  }
  for (int32_t i = 0; i < nodes_count_; ++i) {
    const knhx1_t* node = node_by_id(i);
    id_by_name_.emplace(std::pair<std::string, int32_t>(node->name, i));
  }
  upper_nodes_.resize(nodes_count_);
  fill_upper_nodes(root_);
  fill_tags(buffer.str());
}

void TreeIndex::fill_tags(const std::string& newick_string) {
  tags_.resize(nodes_count_);
  auto parts_by_comma = split(newick_string, ',');
  std::vector<std::string> parts;
  for (auto& part : parts_by_comma) {
    for (auto& part_by_bracket : split(part, ')')) {
      parts.push_back(part_by_bracket);
    }
  }
  for (int32_t id = 0; id < nodes_count_; ++id) {
    auto& part = parts[id];
    if (part.find('[') == std::string::npos) {
      continue;
    }
    std::string tag_string = part.substr(part.find('[') + 1, part.length() - part.find('[') - 2);
    auto tags = split(tag_string, ':');
    for (auto& tag : tags) {
      auto pair = split(tag, '=');
      if (pair.size() != 2) {
        continue;
      }
      tags_[id].emplace(std::pair<std::string, std::string>(pair[0], pair[1]));
    }
  }
  joined_tags_.resize(nodes_count_);
  for (int32_t id = 0; id < nodes_count_; ++id) {
    std::string joined_tags;
    auto& tags = tags_[id];
    auto sci_name = tags.find("sci_name");
    if (sci_name != tags.cend()) {
      joined_tags += "\tsn:Z:" + (*sci_name).second;
    }
    auto rank = tags.find("rank");
    if (rank != tags.cend()) {
      joined_tags += "\tra:Z:" + (*rank).second;
    }
    joined_tags_[id] = joined_tags;
  }

//  for (int32_t id = 0; id < nodes_count_; ++id) {
//    std::cerr << id << " ";
//    for (auto& pair : tags_[id]) {
//      std::cerr << pair.first << ":" << pair.second << " ";
//    }
//    std::cerr << std::endl;
//  }
}

void TreeIndex::fill_upper_nodes(const knhx1_t* node) {
  int32_t node_id = id_by_name(node->name);
  // std::cerr << "processing node: " << node_id << " " << node->name << " " << node->n << std::endl;
  for (int32_t i = 0; i < node->n; ++i) {
    auto child = node_by_id(node->child[i]);
    upper_nodes_[node->child[i]] = upper_nodes_[node_id];
    upper_nodes_[node->child[i]].push_back(node_id);
    fill_upper_nodes(child);
  }
}

std::ostream& operator<<(std::ostream& sout, const TreeIndex& tree) {
  kstring_t str;
  for (int32_t i = 0; i < tree.nodes_count(); ++i) {
    const knhx1_t *p = tree.node_by_id(i);
    printf("[%d] %s\t%d\t%d\t%g", i, p->name, p->parent, p->n, p->d);
    for (int32_t j = 0; j < p->n; ++j) {
      printf("\t%d", p->child[j]);
    }
    putchar('\n');
  }
  str.l = str.m = 0; str.s = 0;
  kn_format(tree.first_node(), tree.nodes_count() - 1, &str);
  puts(str.s);
  free(str.s);
  return sout;
}
