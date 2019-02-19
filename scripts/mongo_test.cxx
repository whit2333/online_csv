#include <fmt/core.h>
#include <fmt/ostream.h>
R__LOAD_LIBRARY(libfmt.so)

#ifdef __cpp_lib_filesystem
#include <filesystem>
namespace fs = std::filesystem;
#else
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
#endif

#include "nlohmann/json.hpp"

#include <iostream>

R__LOAD_LIBRARY(libbsoncxx.so)
R__LOAD_LIBRARY(libmongocxx.so)
#include <bsoncxx/builder/stream/document.hpp>
#include <bsoncxx/json.hpp>

#include <mongocxx/client.hpp>
#include <mongocxx/instance.hpp>

// int main(int, char**) {
void mongo_test() {
  using nlohmann::json;
  json j;
  {
    std::ifstream json_input_file("db2/run_list_coin.json");
    json_input_file >> j;
  }
  mongocxx::instance inst{};
  mongocxx::client   conn{mongocxx::uri{"mongodb://cdaql1.jlab.org:27017"}};

  bsoncxx::builder::stream::document document{};

  auto collection = conn["testdb"]["run_list"];
  //document << "hello"
  //         << "world";
  //collection.insert_one(document.view());
  //collection.insert_one(bsoncxx::from_json(j.dump()));

  for (json::iterator it = j.begin(); it != j.end(); ++it) {
    json j_entry;
    j_entry["run_number"] = it.key();
    for (json::iterator it2 = it.value().begin(); it2 != it.value().end(); ++it2) {
      j_entry[it2.key()] = it2.value();
      //std::cout << it2.key() << " : " << it2.value() << "\n";
    }
    collection.insert_one(bsoncxx::from_json(j_entry.dump()));
  }

  //auto cursor = collection.find({});
  //for (auto&& doc : cursor) {
  //  std::cout << bsoncxx::to_json(doc) << std::endl;
  //}
}

