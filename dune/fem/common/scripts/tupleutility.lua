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
-- Date: 2013-03-09
----------------------------------------------------------------------

dofile( "common.lua" )

function tuple_push_back ( i )
  body = "template< class T" .. i .. ( ( (i == 1) ) and "" or ( ", " .. general_list( i-1, "class T" ) ) ) .. " >\n"
  body = body .. "inline " .. tuple( i ) .. " tuple_push_back ( const " .. tuple( i-1 ) .. " &t, T" .. i .. " t" .. i .. " )" .. "\n"
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
  body = "template< " .. general_list( i, "class T" ) .. " >\n"
  body = body .."inline " .. tuple( i ) .. " tuple_push_front ( const " .. tuple( i-1, 2 ) .. " &t, T1 t1 )" .. "\n"
  body = body .. "{" .. "\n"
  if i == 1 then
    body = body .. "  return " .. tuple( i ) .. "( t1 );" .. "\n"
  else
    body = body .. "  return " .. tuple( i ) .. "( t1, " .. tuple_get_elements( i-1, 1 ) .. " );" .. "\n"
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
