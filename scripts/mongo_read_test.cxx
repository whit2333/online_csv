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
void mongo_read_test() {
  using nlohmann::json;
  mongocxx::instance inst{};
  mongocxx::client   conn{mongocxx::uri{"mongodb://cdaql1.jlab.org:27017"}};

  bsoncxx::builder::stream::document document{};

  auto collection = conn["testdb"]["run_list"];
  
  json j = {{"run_number", "7233"}};
  
  // bsoncxx::stdx::optional<bsoncxx::document::value> maybe_result
  auto maybe_result = collection.find_one(bsoncxx::from_json(j.dump()));
  if(maybe_result) {
    // Do something with *maybe_result;
    std::cout << bsoncxx::to_json(*maybe_result) << "\n";
  }
}

