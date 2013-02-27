----------------------------------------------------------------------
-- File: tupleutility.lua
--
-- Description:
--   This script writes out the definition of Dune::tuple related
--   functions:
--     tuple_push_back
--     tuple_push_front
--     tuple_pop_back
--     tuple_pop_front
--
-- Author: Christoph Gersbacher
-- Date: 2013-02-24
----------------------------------------------------------------------

function template_list ( length, begin )
  begin = ( begin ) and begin or 1 -- use default value begin = 1

  if length == 0 then
    return "tuple<>"
  end

  line = "template< class T" .. begin
  for j = (begin+1), (begin+length-1) do
    line = line .. ", class T" .. j
  end
  line = line .. " >"
  return line
end

function tuple ( length, begin )
  begin = ( begin ) and begin or 1 -- use default value begin = 1

  if length == 0 then
    return "tuple<>"
  end

  line = "tuple<"
  last = (begin+length-1)
  for j = begin, last do
    line = line .. " T" .. j
    line = line .. ( ( j == last ) and " " or "," )
  end
  line = line .. ">"
  return line
end

function tuple_get_elements ( length, begin )
  begin = ( begin ) and begin or 1 -- use default value begin = 1

  if length == 0 then
    return ""
  end

  line = ""
  last = (begin+length-1)
  for j = begin, last do
    line = line .. "get< " .. j-1 .. " >( t )"
    line = line .. ( ( j == last ) and "" or ", " )
  end
  return line
end

function tuple_push_back ( i )
  body = template_list( i ).. "\n"
  body = body .. "inline " .. tuple( i ) .. " tuple_push_back ( const " .. tuple( i-1 ) .. " &t, const T" .. i .. " &t" .. i .. " )" .. "\n"
  body = body .. "{" .. "\n"
  if i == 1 then
    body = body .. "  return " .. tuple( i ) .. "( t" .. i .. " );" .. "\n"
  else
    body = body .."  return " .. tuple( i ) .. "( " .. tuple_get_elements( i-1 ) .. ", t" .. i .. " );" .. "\n"
  end
  body = body .. "}"
  return body
end

function tuple_push_front ( i )
  body = template_list( i ) .. "\n"
  body = body .."inline " .. tuple( i ) .. " tuple_push_front ( const " .. tuple( i-1, 2 ) .. " &t, const T1 &t1 )" .. "\n"
  body = body .. "{" .. "\n"
  if i == 1 then
    body = body .. "  return " .. tuple( i ) .. "( t" .. i .. " );" .. "\n"
  else
    body = body .. "  return " .. tuple( i ) .. "( t1, " .. tuple_get_elements( i-1, 2 ) .. " );" .. "\n"
  end
  body = body .. "}"
  return body
end

function tuple_pop_back ( i )
  body = template_list( i ) .. "\n"
  body = body .. "inline " .. tuple( i-1 ) .. " tuple_pop_back ( const " .. tuple( i ) .. " &t )" .. "\n"
  body = body .. "{" .. "\n"
  if i == 1 then
    body = body .. "  return " .. tuple( i-1 ) .. "();" .. "\n"
  else
    body = body .. "  return " .. tuple( i-1 ) .. "( " .. tuple_get_elements( i-1 ) .. " );" .. "\n"
  end
  body = body .. "}"
  return body
end

function tuple_pop_front ( i )
  body = template_list( i ) .. "\n"
  body = body .. "inline " .. tuple( i-1, 2 ) .. " tuple_pop_front ( const " .. tuple( i ) .. " &t )" .. "\n"
  body = body .. "{" .. "\n"
  if i == 1 then
    body = body .. "  return " .. tuple( i-1, 2 ) .. "();" .. "\n"
  else
    body = body .. "  return " .. tuple( i-1, 2 ) .. "( " .. tuple_get_elements( i-1, 2 ) .. " );" .. "\n"
  end
  body = body .. "}"
  return body
end
