file(REMOVE_RECURSE
  "python/nanogui.pdb"
  "python/nanogui.cpython-36m-darwin.so"
)

# Per-language clean rules from dependency scanning.
foreach(lang C CXX)
  include(CMakeFiles/nanogui-python.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
