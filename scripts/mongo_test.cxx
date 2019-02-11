#include <iostream>

R__LOAD_LIBRARY(libbsoncxx.so)
R__LOAD_LIBRARY(libmongocxx.so)
#include <bsoncxx/builder/stream/document.hpp>
#include <bsoncxx/json.hpp>

#include <mongocxx/client.hpp>
#include <mongocxx/instance.hpp>

//int main(int, char**) {
void mongo_test(){
    mongocxx::instance inst{};
    mongocxx::client conn{mongocxx::uri{"mongodb://cdaql1.jlab.org:27017"}};

    bsoncxx::builder::stream::document document{};

    auto collection = conn["testdb"]["testcollection"];
    document << "hello" << "world";

    collection.insert_one(document.view());
    auto cursor = collection.find({});

    for (auto&& doc : cursor) {
        std::cout << bsoncxx::to_json(doc) << std::endl;
    }
}

