--[[ This script creates an AquaME model, i.e., a TerraME model based on the AQUA algorithm described in Rinaldi et al., 2007.
  AquaME is tried with a scenario SimplePlain that corresponds to the "synthetic scenario" mentioned in that paper.

  Rinaldi, P.R., Dalponte, D.D., Vénere, M.J., Rinaldi, A., 2007. Cellular automata algorithm for simulation of surface flows in large plains.
  Simul. Model. Pract. Theory 15 (3), 315–327.

  "--**" and "--[[**" inhibit code useful for testing and time measuring.
--]]

print(sessionInfo().time, ": Session starts:")

AquaME = Model {
  -- Default name for instances of the model.
  name = "Empty plain",
  -- Default maximum duration, time steps.
  duration = Choice {min = 0.001}, -- hours
  -- Default step duration, in hours.
  stepHours = Choice {min = 0.001}, -- hours per time step.
  -- Default and constraints for sizes of the cell grid:
  xCellSize = Choice {min = 0.001}, -- cells
  yCellSize = Choice {min = 0.001}, -- cells
  -- Default and constraints for actual lenght of the area modeled, in meters:
  xAreaSize = Choice {min = 0.001}, -- meters
  yAreaSize = Choice {min = 0.001}, -- meters
  -- Default and constraints for initial relaxation parameter for all cells, representing the flow resistance.
  alpha = Choice {min = 0, max = 1, default = 0.988},
  -- Default and constraints for alpha0, a constant for calculating alpha for rivers.
  alpha0 = Choice {min = 0, max = 1, default = 0.012},
  -- Default and constraints for "n", a constant that is exponent for calculation of alpha for rives.
  nExp = Choice {min = 0, max = 3, default = 3},
  -- Default and constraints for coefficient representing soil saturation characteristics.
  beta = Choice {min = 0, max = 1, default = 0.9999},
  -- Default and minimum bias infiltration for cells.
  infiltrationBias = Choice {min = 0}, -- meters/hours

  -- Load scenario terrain data.
  loadTerrain = function (_,_) end,

  -- Load scenario precipitation data.
  loadPrecipitation = function (_,_,_) return 0 end,

  -- Return optional parms specific for instance.
  otherParms = function (_) return end,

  -- Initialize model:
  init = function(model)

    -- Maximum duration, time steps.
    model.finalTime = math.ceil(model.duration / model.stepHours) + 1

    model.xdim = math.ceil(model.xAreaSize / model.xCellSize) -- Number or cells in the x dimension.
    model.ydim = math.ceil(model.yAreaSize / model.yCellSize) -- Number or cells in the y dimension.

    -- Sequence of the 9 partitions for runoff calculation (Rinaldi at al, Figure 5).
    model.sequence = {{2,2},{1,1},{2,1},
                      {1,2},{2,3},{1,3},
                      {3,2},{3,1},{3,3}}

    -- No outflow at the beginning of the simulation.
    model.outflow = 0

    model.cell = Cell {
      init = function(cell)
        -- Set defaults for cell using nil for lower memory consumption:
        -- cell.open               -- Cell has no open border.
        -- cell.isRiver            -- Cell is not on a river.
        -- cell.infiltration       -- No initial infiltration.
        -- cell.infiltrationBias   -- Typical infiltration bias for terrain from model.
        -- cell.alpha = model.alpha  -- Typical alpha for terrain from model.

        -- Set defaults for cell using atibuttes because they are usually updated in each time step:
        cell.height = 0 -- All plain.
        cell.water = 0  -- No initial water.

        -- Fix cells' values from default ones according to scenario:
        model:loadTerrain(cell)

      end
    }
    --** print(sessionInfo().time, ": Creates cellular space.")
    model.cs = CellularSpace {xdim = model.xdim, ydim = model.ydim, instance = model.cell}

    --** print(sessionInfo().time, ": Creates cellular space neighborhood with Moore strategy (3x3 grid).")
    model.cs:createNeighborhood {}

    --** print(sessionInfo().time, ": Creates 9 trajectories representing the 9 partitions with their grids.")
    model.partition = {};

    for t = 1, #model.sequence do

      local x = model.sequence[t][1] -- The x reference for the center cells of this partition.
      local y = model.sequence[t][2] -- The y reference for the center cells of this partition.

      model.partition[t] = Trajectory {
        -- Each partition contains cells defined by (x,y) in such a way that they are separated by 2 cells above and below, 
        -- so there is no neighbor at the same partition.
        target = model.cs, select = function(cell) return (1 + cell.x % 3) == x and (1 + cell.y % 3) == y end
      }
    end

    --** print(sessionInfo().time, ": Creates grids around each cell.")
    model.grid = {}
    forEachCell(model.cs, function(cell)
      -- Creates a grid around the cell:
      model.grid[cell] = {cell}
      -- Add cell's neighbors to the grid:
      forEachNeighbor(cell, function(neighbor) table.insert(model.grid[cell], neighbor) end)
      -- Order the cells index in the grid according to the terrain height, in ascending order (Rinaldi et al. 2007, equation 1).
      table.sort(model.grid[cell], function(lower, higher) return lower.height < higher.height end)
    end)

    --** print(sessionInfo().time, ": Define chart.")

    -- [[** To turning chart off.
    model.chart = Chart {
      title = model.name,
      target = model, --** Alternative: model.cs:get(2, model.ydim)
      select = "outflow", --** Alternatives: water, infiltration, alpha, height.
      color = "blue"
    }
    --]]

    --** print(sessionInfo().time, ": Define timed events.")
    model.timer = Timer {

      Event {priority = 1, action = function(event)

        --** print(sessionInfo().time, ": Calculates precipitation, infiltration and alpha for rivers.")

        forEachCell(model.cs, function(cell)

          -- Adds precipitation to cell water, if there is any defined in model.
          cell.water = cell.water + model:loadPrecipitation(cell, event:getTime())

          --** print(sessionInfo().time, ": Calculates cell infiltration.")

          -- Using nil for defaults for lower memory consumption.
          local infiltration = cell.infiltration or 0
          local infiltrationBias = cell.infiltrationBias or model.infiltrationBias
          local beta = cell.beta or model.beta

          if cell.water < beta * infiltration + infiltrationBias
            then
              cell.infiltration = cell.water
              cell.water = 0
            else
              cell.infiltration = beta * infiltration + infiltrationBias
              cell.water = cell.water - cell.infiltration
          end

          -- Alpha for rivers are calculated from its current water level (Rinaldi et al. 2007, equation 8)
          if cell.isRiver then
            local alpha = 1 - model.alpha0 * cell.water ^ model.nExp
            -- Alpha can't be lower than zero.
            if alpha < 0 then alpha = 0 end
            --[[** Alpha for rivers shouldn't be higher than alpha for terrain, but authors oriented this way:
            if alpha > model.alpha then alpha = model.alpha end
            --]]
            -- Alpha is recalculated and saved for rivers only: for other cells, alpha is not saved and always assumed to be model.alpha.
            cell.alpha = alpha
          end
        end)
      end},

      Event {priority = 2, action = function(_)

        --** print(sessionInfo().time, ": Calculates runoff.")

        -- For all partitions ordered by model.sequence (Rinaldi et al. 2007, figure 5).
        for t = 1, #model.sequence do

          forEachCell(model.partition[t], function(cell)

            -- Local variable for improved readability and performance:
            local grid = model.grid[cell] -- Obtain the gridaroung the cell previously ordered by height.
            local N = #grid -- The size of the grid ordered by height. Mostly 3x3 = 9, but could be 4 at corners and 6 at borders.
            local wOld = {} -- Current water level of the ith cell of the grid ordered by height.
            local h = {}    -- Cell height of the ith cell of the ordered grid ordered by height.
            local W = 0     -- Current water volume contained in the grid divided by the unit-cell area.

            --For each ith cell in the ordered grid:
            for i = 1, N do
              wOld[i] = grid[i].water; h[i] = grid[i].height
              -- W is the current water volume contained in the grid divided by the unit-cell area (Rinaldi et al. 2007, equation 3).
              W = W + wOld[i]
            end

            -- If there is no water in the grid, there is nothing to calculate.
            if W == 0 then return end

            -- Obtain k, the number of cells that remains wet after the water drains down (Rinaldi et al. 2007, equation 4):
            local k = N
            repeat
              local s = 0; for i = 1, k do s = s + h[k] - h[i] end
              k = k - 1
            until W >= s or k == 0
            k = k + 1

            -- Calculate the equilibrium surface water height, H, as the height the water would reach (Rinaldi et al. 2007, equation 2):
            local H = W/k; for i = 1, k do H = H + h[i]/k end

            -- Calculate the new water level of ith cell of the grid as a linear combination of wOld[i] and the drained water level (Rinaldi et al. 2007, equation 6):
            -- The drained water level wDrain[i] is zero for cells below k and H - h[i] for cells above (Rinaldi et al. 2007, equation 5).
            -- The value of alpha corresponds to the center cell of the grid.
            local alpha = cell.alpha or model.alpha -- Using nil as default for lower memory consumption.
            for i = 1, k      do grid[i].water = alpha * wOld[i] + (1 - alpha) * (H - h[i]) end
            for i = k + 1, N  do grid[i].water = alpha * wOld[i] end

            --[[** Checks conservation of water at the grid:
            local Wnew = 0; for i = 1, k do Wnew = Wnew + H - h[i] end
            if W > 0 and (W - Wnew)/W > 10^(-8) then
              print("W changed at: (",cell.x,",",cell.y,") in:", 100*(Wnew-W)/W, "%.")
            end
            --**]]

          end)
        end
      end},

      Event {priority = 3, action = function(_)

        --** print(sessionInfo().time, ": Calculates outflow.")

        local outflow = 0

        forEachCell(model.cs, function(cell)

          -- If cell at an open border,let water ouflow from cell .
          if cell.open then outflow = outflow + cell.water; cell.water = 0 end
        
        end)

        -- Calculates outflow in m3/sec:
        model.outflow = outflow
          * model.xCellSize  -- meters per cell, x dimension
          * model.yCellSize  -- meters per cell, y dimension
          / model.stepHours  -- hours per time step
          / 3600             -- seconds per hour

      end},

      -- Updates chart.
      Event {priority = 4, action = model.chart},

      --[[**
      Event {priority = 5, period = 1000, action = function(event)
        print(sessionInfo().time, ": End of time step ", event:getTime(), ".")
      end
      --]]
    }
  end
}
--** print(sessionInfo().time, ": AquaME base model created.")

-- Creates a AquaME scenario redefining model parameters and functions.
simplePlain = AquaME {
  name = "Simple plain",
  duration = 40,           -- Duration of the simulation in hours.
  stepHours = 1/1000,      -- Hours per time step.
  xCellSize = 80,          -- x cell size in meters. xdim adjusted to ceil.
  yCellSize = 80,          -- y cell size in meters. ydim adjusted to ceil.
  xAreaSize = 11000,       -- Size of the x dimension in meters.
  yAreaSize = 8800,        -- Size of the y dimension in meters.
  alpha = 0.988,           -- Flow resistance for terrain.
  alpha0 = 0.012,          -- Constant for calculating alpha for rivers.

  -- Add parameters for functions specific for scenario.
  otherParms = function (_) return {
    plainSlope = 0.0085,   -- Plain slope at the ydim only, in meters per meter.
    channelSlope = 0.0025, -- Channel slope at the xdim only, in meters per meter.
    channelWidth = 5,      -- Channel width, in meters. Adjusted to grids (ceil).
                           -- If absent or zero, no channel and opens the lower border of plain for outflow.
    channelDepth = 2.5,    -- Additional depth at the highest border of the channel, in meters.
    rainBegin = 0,         -- Instant when when rain begins after the beginning of the simulatio, in hours.
    rainDuration = 10,      -- Rain duration, in hours.
    rainRate = 0.015,      -- Rate of rain, in mm per hour.
    rainAll = false        -- True if it rains also on the channel.
    }
  end,

  -- Loads height and defines rivers for scenario.
  loadTerrain = function(model, cell)

    -- Other parms specific for scenario. Use of local parms for improved performance.
    local parm = model:otherParms()

    -- Defines cell.heigh and channel width (Rinaldi, 2007, figure 14).

    -- Calculates channel width in number of cells.
    local channelGridWidth = 0
    if not parm.channelWidth then
      -- If no channel is defined, than lower border of the plain is all open for outflow.
      cell.open = (cell.x == 0)
    elseif parm.channelWidth == 0 then
      -- If channel defined as zero, then assumes 1 grid.
        channelGridWidth = 1
    elseif parm.channelWidth ~= 0 then
      -- If channel defined and not zero, assumes an integer number of grids (ceil).
      channelGridWidth = math.ceil(parm.channelWidth / model.yCellSize)
    end

    -- If cell is on the channel:
    if cell.y >= model.ydim - channelGridWidth then
      -- Channel behave as a river.
      cell.isRiver = true
      -- Height of cell in the channel calculated from its slope on xdim and doesn't depend on cell.y
      cell.height = cell.x * parm.channelSlope * model.xCellSize
      -- In this scenario, the lowest end of the channel is an open border.
      cell.open = (cell.x == 0)

    -- If cell is not on the channel, cell is on the plain:
    else cell.height =
      -- Plain reachs its lowest level at the border of the channel and doesn't depende on cell.x:
      (model.ydim - 1 - channelGridWidth - cell.y) * parm.plainSlope * model.yCellSize +
      -- The lowest plain border is above the highest channel border...
      (model.xdim - 1) * parm.channelSlope * model.xCellSize +
      -- ... and an additional depth at the highest border of the channel may be present.
      (parm.channelDepth or 0)
    end
  end,

  -- Loads precipitation for the scenario.
  loadPrecipitation = function (model, cell, t)

    -- Other parms specific for scenario. Use of local parms for improved performance.
    local parm = model:otherParms()

    -- Calculate channel width in number of cells.
    local channelGridWidth = 0
    if not parm.channelWidth then
      -- If no channel is defined, than lower border of the plain is all open for outflow.
      cell.open = (cell.x == 0)
    elseif parm.channelWidth == 0 then
      -- If channel defined as zero, then assumes 1 grid.
      channelGridWidth = 1
    elseif parm.channelWidth ~= 0 then
      -- If channel defined and not zero, assumes an integer number of cells (ceil).
      channelGridWidth = math.ceil(parm.channelWidth / model.yCellSize)
    end

    -- Calculate rain begin and end in time steps.
    local rainBegin = parm.rainBegin // model.stepHours + 1
    local rainEnd = rainBegin - 1 + parm.rainDuration // model.stepHours

    -- Calculate rain rate in mm per time step.
    local rainRate = parm.rainRate * model.stepHours

    if t >= rainBegin and t <= rainEnd and
      -- Depending on a scenario option, it rains all over the scenario, ...
      (parm.rainAll or
      -- ... or only on the plain cells
      cell.y <= model.ydim - 1 - channelGridWidth)
      -- Fixed rain rate during the period for all cells in the region:
      then return rainRate
      -- No rain out of the time interval:
      else return 0
    end
  end,
}

--** print(sessionInfo().time, ": SimplePlain instance created.")
print(sessionInfo().time, ": Run scenario.")

simplePlain:run()

print(sessionInfo().time, ": Session ends.")