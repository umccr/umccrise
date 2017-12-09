function count_true(...)
  local t = {...}
  local count = 0
  for i, v in pairs(t) do
    if v == true then
      count = count + 1
    end
  end
  return count
end