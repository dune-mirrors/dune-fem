----------------------------------------------------------------------
-- File: streams.lua
--
-- Description:
--   Please doc me.
--
-- Author: Christoph Gersbacher
-- Date: 2013-03-09
----------------------------------------------------------------------

dofile( "common.lua" )

function template_stream_list ( length, begin )
  if length == 0 then
    return "template< " .. begin .. " >"
  else 
    prefix = "class T"
    return "template< " ..begin .. ", " .. general_list( length, prefix ) .. " >"
  end
end

function tuple_to_instream ( i )
  body = template_stream_list( i, "class StreamTraits" ) .. "\n"
  body = body .. "inline InStreamInterface< StreamTraits > &\n"
  body = body .. "  operator>> ( InStreamInterface< StreamTraits > &in, " .. tuple( i ) .. " &tuple )\n"
  if i == 0 then
    body = body .. "{\n  return in;\n}"
    return body
  end
  body = body .. "{\n"
  body = body .. "  typedef " .. tuple( i ) .. " Tuple;\n"
  body = body .. "  return TupleToInStream< InStreamInterface< StreamTraits > >::template apply< Tuple >( in, tuple );\n"
  body = body .. "}"
  return body
end

function tuple_to_outstream ( i )
  body = template_stream_list( i, "class StreamTraits" ) .. "\n"
  body = body .. "inline OutStreamInterface< StreamTraits > &\n"
  body = body .. "  operator<< ( OutStreamInterface< StreamTraits > &out, const " .. tuple( i ) .. " &tuple )\n"
  if i == 0 then
    body = body .. "{\n  return out;\n}"
    return body
  end
  body = body .. "{\n"
  body = body .. "  typedef " .. tuple( i ) .. " Tuple;\n"
  body = body .. "  return TupleToOutStream< OutStreamInterface< StreamTraits > >::template apply< Tuple >( out, tuple );\n"
  body = body .. "}"
  return body
end

local tuple_length = 9

print "// Dune::tuple to InStream" 
print "// -----------------------"

print ""

for i = 0, tuple_length do
  print( tuple_to_instream( i ) )
  print ""
end

print ""
print ""

print "// Dune::tuple to OutStream" 
print "// ------------------------"

print ""

for i = 0, tuple_length do
  print( tuple_to_outstream( i ) )
  print ""
end
