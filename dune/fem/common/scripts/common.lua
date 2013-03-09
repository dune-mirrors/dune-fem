----------------------------------------------------------------------
-- File: common.lua
--
-- Description:
--   Please doc me.
--
-- Author: Christoph Gersbacher
-- Date: 2013-03-09
----------------------------------------------------------------------

function general_list ( length, prefix, begin )
  begin = ( begin ) and begin or 1      -- use default value begin = 1
  prefix = ( prefix ) and prefix or "T"  -- use default value ""

  if length == 0 then
    return ""
  end

  line = prefix .. begin
  for j = (begin+1), (begin+length-1) do
    line = line .. ", " .. prefix .. j
  end
  return line
end

function template_list ( length, begin )
  prefix = "class T"
  if length == 0 then
    return "template<>"
  else 
    return "template< " .. general_list( length, prefix, begin ) .. " >"
  end
end

function tuple ( length, begin )
  if length == 0 then
    return "Dune::tuple<>"
  else 
    return "Dune::tuple< " .. general_list( length, "T", begin ) .. " >"
  end
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
