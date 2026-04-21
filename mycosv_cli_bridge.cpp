#define main mycosv_cli_main
#include "main.cpp"
#undef main

extern "C" {
int run_mycosv_cli(int argc, const char** argv) {
    return mycosv_cli_main(argc, const_cast<char**>(argv));
}
}
