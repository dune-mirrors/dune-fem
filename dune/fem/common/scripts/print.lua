dofile( "tupleutility.lua" )

local tuple_length = 9

print "// tuple_push_back"
print "// ---------------"

print ""

for i = 1, tuple_length do
  body = tuple_push_back( i )
  print( body )
  print ""
end

print ""
print ""

print "// tuple_push_front"
print "// ----------------"

print ""

for i = 1, tuple_length do
  body = tuple_push_front( i )
  print( body )
  print ""
end

print ""
print ""

print "// tuple_pop_back"
print "// ----------------"

print ""

for i = 1, tuple_length do
  body = tuple_pop_back( i )
  print( body )
  print ""
end

print ""
print ""

print "// tuple_pop_front"
print "// ----------------"

print ""

for i = 1, tuple_length do
  body = tuple_pop_front( i )
  print( body )
  print ""
end
