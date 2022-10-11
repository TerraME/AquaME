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
  -- Neighborhood strategy.
  strategy = Choice {"moore", "vonneumann"},
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
  loadTerrain = function () end,

  -- Load scenario precipitation data.
  loadPrecipitation = function () return 0 end,

  -- Return optional parms specific for scenario instance.
  loadScenarioParms = function () end,

  -- Initialize model:
  init = function(model)

    -- Loads into model any addtional parm needed for scenario instance:
    model:loadScenarioParms()

    -- Calculates the maximum duration in time steps:
    model.finalTime = math.ceil(model.duration / model.stepHours) + 1

    -- Calculates the number or cells in the x and y dimensions:
    model.xdim = math.ceil(model.xAreaSize / model.xCellSize)
    model.ydim = math.ceil(model.yAreaSize / model.yCellSize)

    -- Define the sequence of partitions for runoff calculation according to the strategy and to the center cell x and y coordinates.
    if model.strategy == "vonneumann"
      then model.sequence = { {1,2,3,4,5},
                              {4,5,1,2,3},
                              {2,3,4,5,1},
                              {5,1,2,3,4},
                              {3,4,5,1,2}
                            }
      else model.sequence = { {2,4,6}, -- (Rinaldi at al, Figure 5)
                              {3,1,5},
                              {8,7,9}
                            }
    end

    -- No outflow at the beginning of the simulation.
    model.outflow = {}

    model.cell = Cell {
      init = function(cell)
        -- Many defaults use nil for lower memory consumption:
        -- cell.open                -- Cell has no open border.
        -- cell.isRiver             -- Cell is not on a river.
        -- cell.infiltration        -- No initial infiltration.
        -- cell.infiltrationBias    -- Typical infiltration bias for terrain comes from model.
        -- cell.alpha               -- Typical alpha for terrain comes from model.

        -- Set defaults for cell using atibuttes because they are usually updated in each time step:
        cell.height = 0 -- All plain.
        cell.water = 0  -- No initial water.

        -- Fix cells' height and alpha from default ones according to scenario:
        model:loadTerrain(cell)

      end
    }
    --** print(sessionInfo().time, ": Creates cellular space.")
    model.cs = CellularSpace {xdim = model.xdim, ydim = model.ydim, instance = model.cell}

    --** print(sessionInfo().time, ": Creates cellular space neighborhood with model strategy.")
    model.cs:createNeighborhood {strategy = model.strategy}

    --** print(sessionInfo().time, ": Creates a grid around each cell.")
    model.grid = {}
    forEachCell(model.cs, function(cell)

      -- Calculates the index (i.e., the order of evaluation) of the partition that the grid around the cell belongs to:
      local partition = model.sequence[1 + cell.x % #model.sequence][1 + cell.y % #model.sequence]

      -- Creates the partition, if it hasn't been created yet:
      model.grid[partition] = model.grid[partition] or {}

      -- Initiates the grid around the cell with the cell itself, on the partition the grid belongs to:
      model.grid[partition][cell] = {cell}

      -- Adds cell's neighbors to this grid:
      forEachNeighbor(cell, function(neighbor) table.insert(model.grid[partition][cell], neighbor) end)

      -- Orders the cells in the grid according to the terrain height, in ascending order (Rinaldi et al. 2007, equation 1).
      table.sort(model.grid[partition][cell], function(lower, higher) return lower.height < higher.height end)

    end)

    model.timer = Timer {

        -- Calculates outflow, precipitation, infiltration and alpha for rivers.
        Event {priority = 1, action = function(event)

        -- Initialize event outflow with zero.
        local outflow = 0

        -- Update water and outflow of each cell for the event.
        forEachCell(model.cs, function(cell)

          -- If cell is at an open border, let water ouflow from cell.
          if cell.open then outflow = outflow + cell.water; cell.water = 0 end

          -- Adds precipitation to cell water.
          cell.water = cell.water + model:loadPrecipitation(cell, event:getTime())

          -- Replace nil by default values:
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
            --[[** Alpha for rivers shouldn't be higher than alpha for terrain, but Rinaldi et al. (2007) oriented this way:
            if alpha > model.alpha then alpha = model.alpha end
            --]]
            cell.alpha = alpha
          end
        end)

        -- Calculates and save outflow in m3/sec for the time step:
        model.outflow[event:getTime()-1] = outflow
          * model.xCellSize  -- meters per cell, x dimension
          * model.yCellSize  -- meters per cell, y dimension
          / model.stepHours  -- hours per time step
          / 3600             -- seconds per hour

        -- Calculates runoff for each partition (Rinaldi et al. 2007, figure 5):
        for partition = 1, #model.grid do

          -- For each grid of the partition, obtains the center cell and its grid:
          for cell, grid in pairs(model.grid[partition]) do
            --[[ The steps in this loop were carefully crafted in order to enhance performance,
              as this is the most repeated loop in the algorithm.
            --]]

            local N = #grid                         -- The size of the grid. Depends on the strategy and if the cell is at a corner or border.
            local h = {}; h[N] = grid[N].height     -- Will store the height of each cell.
            local s = 0                             -- Will store the sum of h[i] from h[1] to h[N-1]

            local W = grid[N].water                 -- Will store the current water volume contained in the grid divided by the unit-cell area  (Rinaldi et al. 2007, equation 3).
            local alpha = cell.alpha or model.alpha -- Resistance to runoff. Using nil as default for lower memory consumption.
            grid[N].water = alpha * grid[N].water   -- Updates the next water level. Adjusted below for cells that receive water from higher cells.

            for i = 1, N-1 do
              h[i] = grid[i].height
              s = s + h[i]
              W = W + grid[i].water
              grid[i].water = alpha * grid[i].water
            end

            -- If there is no water in the grid, there is nothing else to calculate,
            -- Testing this helps performance when there is a large number of grids with no water every time step.
            if W == 0 then return end

            -- Obtain k, the number of cells that remains wet after the water drains down (Rinaldi et al. 2007, equation 4):
            -- and adjusts "s", the sum of h[i] until h[k-1].
            local k = N
            while W + s < (k-1)*h[k] and k >= 2 do 
              k = k - 1
              s = s - h[k] 
            end

            -- Calculate H, the equilibrium surface height the water would reach if enough time was provided (Rinaldi et al. 2007, equation 2):
            local H = (W + s + h[k])/k

            -- The added drained water level wDrain[i] is H - h[i] for cells below k (Rinaldi et al. 2007, equation 5).
            -- The new water level of ith cell of the grid as a linear combination of current water and drained water levels (Rinaldi et al. 2007, equation 6):
            for i = 1, k do grid[i].water = grid[i].water + (1 - alpha) * (H - h[i]) end

          end
        end
      end},

      --[[**
      Event {priority = 5, period = 5000, action = function(event)
        print(sessionInfo().time, ": End of time step ", event:getTime(), ".")
      end},
      --]]

      -- Shows chart:
      Event {priority = 6, period = model.finalTime, action = function(_)
        Chart{target = DataFrame{outflow = model.outflow}, select = "outflow", color = "blue" }
      end},
    }
  end
}
--** print(sessionInfo().time, ": AquaME base model created.")

-- Creates a AquaME scenario redefining model parameters and functions.
simplePlain = AquaME {
  strategy = "vonneumann",
  name = "Simple plain",
  duration = 40,           -- Duration of the simulation in hours.
  stepHours = 1/1000,      -- Hours per time step.
  xCellSize = 80,          -- x cell size in meters. xdim adjusted to ceil.
  yCellSize = 80,          -- y cell size in meters. ydim adjusted to ceil.
  xAreaSize = 11000,       -- Size of the x dimension in meters.
  yAreaSize = 8800,        -- Size of the y dimension in meters.
  alpha = 0.988,           -- Flow resistance for terrain.
  alpha0 = 0.012,          -- Constant for calculating alpha for rivers.

  -- Adds specific parameters to the scenario.
  loadScenarioParms = function (model)

    model.plainSlope = 0.0085   -- Plain slope at the ydim only, in meters per meter.
    model.channelSlope = 0.0025 -- Channel slope at the xdim only, in meters per meter.
    model.channelWidth = 5      -- Channel width, in meters. It will be adjusted to an integer number of grids (ceil).

    model.channelDepth = 2.5    -- Additional depth at the highest border of the channel, in meters.
    model.rainBegin = 0         -- Instant when when rain begins after the beginning of the simulatio, in hours.
    model.rainDuration = 10     -- Rain duration, in hours.
    model.rainRate = 0.015      -- Rate of rain, in m per hour.
    model.rainAll = false       -- True if it rains also on the channel.

    -- Assumes an integer number of grids for channel greater than zero.
    model.channelGridWidth = math.ceil(model.channelWidth / model.yCellSize)

    -- Calculates rain begin and end in time steps.
    model.stepRainBegin = model.rainBegin // model.stepHours + 1
    model.stepRainEnd = model.stepRainBegin - 1 + model.rainDuration // model.stepHours

    -- Calculates rain rate in mm per time step.
    model.stepRainRate = model.rainRate * model.stepHours

  end,

  -- Loads height and defines rivers for scenario (Rinaldi, 2007, figure 14).
  loadTerrain = function(model, cell)

    -- If cell is on the channel:
    if cell.y >= model.ydim - model.channelGridWidth
      then
        -- Channel behave as a river.
        cell.isRiver = true
        -- Height of cell in the channel calculated from its slope on xdim and doesn't depend on cell.y
        cell.height = cell.x * model.channelSlope * model.xCellSize
        -- The lowest end of the channel is an open border.
        cell.open =(cell.x == 0)
      -- If cell is not on the channel, cell is on the plain:
      else cell.height =
        -- Plain reachs its lowest level at the border of the channel and doesn't depende on cell.x:
        (model.ydim - 1 - model.channelGridWidth - cell.y) * model.plainSlope * model.yCellSize +
        -- The lowest plain border is above the highest channel border...
        (model.xdim - 1) * model.channelSlope * model.xCellSize +
        -- ... and an additional depth at the highest border of the channel may be present.
        (model.channelDepth or 0)
    end

  end,

  -- Loads precipitation for the scenario.
  loadPrecipitation = function (model, cell, t)

    -- Makes it rain during the period set.
    if t <= model.stepRainEnd and t >= model.stepRainBegin and
      -- Depending on a scenario option, it rains all over the scenario or only on the plain.
      (cell.y <= model.ydim - 1 - model.channelGridWidth or model.rainAll)
      -- Fixed rain rate during the period for all cells in the region:
      then return model.stepRainRate
      -- No rain out of the time interval:
      else return 0
    end

  end,

}

-- Runs scenario:
simplePlain:run()

print(sessionInfo().time, ": Session ends.")