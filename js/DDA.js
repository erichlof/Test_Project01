let canvas = document.getElementById("myCanvas");
let context = canvas.getContext("2d");

canvas.width = window.innerWidth;
canvas.height = window.innerHeight;

let debugInfo = document.getElementById("debugInfo");
debugInfo.innerHTML = "This is some debug text";

let cellDisplayScale = 40;
let numOfCellsX = 20;
let numOfCellsY = 15;
let gridSizeX = numOfCellsX * cellDisplayScale;
let gridSizeY = numOfCellsY * cellDisplayScale;
let cellX = 0;
let cellY = 0;
let cellIndex = -1;
let mouseX = -1;
let mouseY = -1;
let mouseCellIndex = -1;
let cellsArray = new Int8Array(numOfCellsX * numOfCellsY);
let leftButtonPressed = false;
let clickStartIndex = -1;
let playerX = Math.random() * gridSizeX * 0.5;
let playerY = Math.random() * gridSizeY * 0.5;
let playerCellX = Math.floor(playerX / cellDisplayScale);
let playerCellY = Math.floor(playerY / cellDisplayScale);
let playerCellIndex = playerCellY * numOfCellsX + playerCellX;
let playerLookDirX = 0;
let playerLookDirY = 0;
let playerLookMagnitude = 0;
let playerLookAngle = 0;
let playerMoveForward = 0; // will be +1 for moving forward, -1 for moving backwards, 0 for no movement
let playerStrafe = 0; // will be +1 for strafing right, -1 for strafing left, 0 for no movement
let testPlayerX = 0;
let testPlayerY = 0;
let testCellX = -1;
let testCellY = -1;
let testCellIndex = -1;
// DDA variables
let rayStartX = 0;
let rayStartY = 0;
let rayDirectionX = 0;
let rayDirectionY = 0;
let rayUnitStepSizeX = 0;
let rayUnitStepSizeY = 0;
let mapCheckX = 0;
let mapCheckY = 0;
let rayLength1D_X = 0;
let rayLength1D_Y = 0;
let stepX = 0;
let stepY = 0;
let tileFound = false;
let rayDistance = 0;
let MAX_RAY_DISTANCE = 1000;
let intersectionX = 0;
let intersectionY = 0;


window.addEventListener("resize", handleWindowResize);
function handleWindowResize()
{
	canvas.width = window.innerWidth;
	canvas.height = window.innerHeight;
}

window.addEventListener("mousedown", handleMouseDown);
function handleMouseDown(event)
{
	event.preventDefault();

	if (event.button == 0)
	{
		leftButtonPressed = true;
		clickStartIndex = mouseCellIndex;

		if (cellX > -1 && cellY > -1)
		{
			if (cellsArray[mouseCellIndex] == 0 && mouseCellIndex != playerCellIndex)
			{
				cellsArray[mouseCellIndex] = 1;
			}
			else
			{
				cellsArray[mouseCellIndex] = 0;
			}
		}
	}
}

window.addEventListener("mouseup", handleMouseUp);
function handleMouseUp(event)
{
	if (event.button == 0)
	{
		leftButtonPressed = false;
		clickStartIndex = -1;
	}
}


window.addEventListener("mousemove", handleMouseMove);
function handleMouseMove(event)
{
	mouseX = event.clientX;
	mouseY = event.clientY;
	cellX = Math.floor(mouseX / cellDisplayScale);
	cellY = Math.floor(mouseY / cellDisplayScale);
	mouseCellIndex = cellY * numOfCellsX + cellX;

	if (cellX < 0 || cellX >= numOfCellsX || cellY < 0 || cellY >= numOfCellsY)
	{
		cellX = cellY = mouseCellIndex = -1;
	}

	if (mouseCellIndex != clickStartIndex && mouseCellIndex != playerCellIndex &&
		cellsArray[mouseCellIndex] == 0 && leftButtonPressed && cellX > -1 && cellY > -1)
	{
		cellsArray[mouseCellIndex] = 1;
	}
}

let keyCode = null;
let KeyA_isPressed = false;
let KeyD_isPressed = false;
let KeyW_isPressed = false;
let KeyS_isPressed = false;
window.addEventListener("keydown", handleKeyDown);
function handleKeyDown(event)
{
	event.preventDefault();

	keyCode = event.code;

	if (keyCode == "KeyA")
	{
		KeyA_isPressed = true;
	}
	if (keyCode == "KeyD")
	{
		KeyD_isPressed = true;
	}
	if (keyCode == "KeyW")
	{
		KeyW_isPressed = true;
	}
	if (keyCode == "KeyS")
	{
		KeyS_isPressed = true;
	}
}

window.addEventListener("keyup", handleKeyUp);
function handleKeyUp(event)
{
	event.preventDefault();

	keyCode = event.code;

	if (keyCode == "KeyA")
	{
		KeyA_isPressed = false;
	}
	if (keyCode == "KeyD")
	{
		KeyD_isPressed = false;
	}
	if (keyCode == "KeyW")
	{
		KeyW_isPressed = false;
	}
	if (keyCode == "KeyS")
	{
		KeyS_isPressed = false;
	}
}


function draw()
{
	// update controls
	playerMoveForward = 0;
	playerStrafe = 0;

	if (KeyW_isPressed && !KeyS_isPressed)
	{
		playerMoveForward = 1;
	}
	if (KeyS_isPressed && !KeyW_isPressed)
	{
		playerMoveForward = -1
	}
	if (KeyD_isPressed && !KeyA_isPressed)
	{
		playerStrafe = 1;
	}
	if (KeyA_isPressed && !KeyD_isPressed)
	{
		playerStrafe = -1;
	}


	playerCellX = Math.floor(playerX / cellDisplayScale);
	playerCellY = Math.floor(playerY / cellDisplayScale);
	playerCellIndex = playerCellY * numOfCellsX + playerCellX;


	// DRAW

	context.clearRect(0, 0, canvas.width, canvas.height);


	for (let row = 0; row < numOfCellsY; row++)
	{
		for (let col = 0; col < numOfCellsX; col++)
		{
			cellIndex = row * numOfCellsX + col;

			context.lineWidth = 1;
			context.strokeStyle = "red";
			context.strokeRect(col * cellDisplayScale, row * cellDisplayScale,
				cellDisplayScale, cellDisplayScale);

			if (cellIndex == mouseCellIndex)
			{
				context.lineWidth = 2;
				context.strokeStyle = "white";
				context.strokeRect(col * cellDisplayScale + 1, row * cellDisplayScale + 1,
					cellDisplayScale - 3, cellDisplayScale - 3);
			}
			if (cellIndex == mouseCellIndex && cellsArray[cellIndex] == 1)
			{
				context.fillStyle = "gray";
				context.fillRect(col * cellDisplayScale + 2, row * cellDisplayScale + 2,
					cellDisplayScale - 5, cellDisplayScale - 5);
			}
			if (cellIndex != mouseCellIndex && cellsArray[cellIndex] == 1)
			{
				context.fillStyle = "gray";
				context.fillRect(col * cellDisplayScale + 1, row * cellDisplayScale + 1,
					cellDisplayScale - 1, cellDisplayScale - 1);
			}
		} // end for (let col = 0; col < numOfCellsX; col++)
	} // end for (let row = 0; row < numOfCellsY; row++)


	// player lookDirection and lookAngle
	playerLookDirX = mouseX - playerX;
	playerLookDirY = mouseY - playerY;
	playerLookMagnitude = Math.sqrt((playerLookDirX * playerLookDirX) + (playerLookDirY * playerLookDirY))
	playerLookDirX /= playerLookMagnitude;
	playerLookDirY /= playerLookMagnitude;
	playerLookAngle = Math.atan2(playerLookDirY, playerLookDirX);

	// DDA algorithm
	rayStartX = playerX / cellDisplayScale;
	rayStartY = playerY / cellDisplayScale;
	mapCheckX = playerCellX;
	mapCheckY = playerCellY;
	rayDirectionX = playerLookDirX;
	rayDirectionY = playerLookDirY;
	rayUnitStepSizeX = Math.sqrt(1 + (rayDirectionY / rayDirectionX) * (rayDirectionY / rayDirectionX));
	rayUnitStepSizeY = Math.sqrt(1 + (rayDirectionX / rayDirectionY) * (rayDirectionX / rayDirectionY));
	if (rayDirectionX < 0)
	{
		stepX = -1;
		rayLength1D_X = (rayStartX - mapCheckX) * rayUnitStepSizeX;
	}
	else
	{
		stepX = 1;
		rayLength1D_X = ((mapCheckX + 1) - rayStartX) * rayUnitStepSizeX;
	}
	if (rayDirectionY < 0)
	{
		stepY = -1;
		rayLength1D_Y = (rayStartY - mapCheckY) * rayUnitStepSizeY;
	}
	else
	{
		stepY = 1;
		rayLength1D_Y = ((mapCheckY + 1) - rayStartY) * rayUnitStepSizeY;
	}

	tileFound = false;
	while (!tileFound)// && rayDistance < MAX_RAY_DISTANCE)
	{
		// Walk
		if (rayLength1D_X < rayLength1D_Y)
		{
			mapCheckX += stepX;
			rayDistance = rayLength1D_X;
			rayLength1D_X += rayUnitStepSizeX;
		}
		else
		{
			mapCheckY += stepY;
			rayDistance = rayLength1D_Y;
			rayLength1D_Y += rayUnitStepSizeY;
		}

		// make sure we're still within grid X and Y boundaries
		if (mapCheckX < 0 || mapCheckX >= numOfCellsX || mapCheckY < 0 || mapCheckY >= numOfCellsY)
		{
			tileFound = true;
			break;
		}

		// Check tile
		if (cellsArray[(mapCheckY * numOfCellsX) + mapCheckX] == 1)
		{
			tileFound = true;
		}
	} // end while (!tileFound)

	if (tileFound)
	{
		intersectionX = playerX + (rayDirectionX * rayDistance * cellDisplayScale);
		intersectionY = playerY + (rayDirectionY * rayDistance * cellDisplayScale);
	}


	if (KeyW_isPressed || KeyS_isPressed)
	{
		testPlayerX = playerX;
		testPlayerY = playerY;
		// move player's X position in the playerLookDir vector, but only use the X direction component
		testPlayerX += playerLookDirX * playerMoveForward;
		//testPlayerY += playerLookDirY * playerMoveForward;

		testCellX = Math.floor(testPlayerX / cellDisplayScale);
		testCellY = Math.floor(testPlayerY / cellDisplayScale);
		testCellIndex = testCellY * numOfCellsX + testCellX;
		if (testCellX > -1 && testCellX < numOfCellsX && testCellY > -1 && testCellY < numOfCellsY &&
			cellsArray[testCellIndex] == 0)
		{ // success with moving in X direction only
			playerX = testPlayerX;
			//playerY = testPlayerY;
		}

		// reset testPlayer X and Y
		testPlayerX = playerX;
		testPlayerY = playerY;
		// move player's Y position in the playerLookDir vector, but only use the Y direction component
		//testPlayerX += playerLookDirX * playerMoveForward;
		testPlayerY += playerLookDirY * playerMoveForward;

		testCellX = Math.floor(testPlayerX / cellDisplayScale);
		testCellY = Math.floor(testPlayerY / cellDisplayScale);
		testCellIndex = testCellY * numOfCellsX + testCellX;
		if (testCellX > -1 && testCellX < numOfCellsX && testCellY > -1 && testCellY < numOfCellsY &&
			cellsArray[testCellIndex] == 0)
		{ // success with moving in Y direction only
			//playerX = testPlayerX;
			playerY = testPlayerY;
		}

	} // end if (KeyW_isPressed || KeyS_isPressed)

	if (KeyA_isPressed || KeyD_isPressed)
	{
		testPlayerX = playerX;
		testPlayerY = playerY;
		// trick to moving perpendicular to playerLookDir vector (sideways strafe)
		// move player's X position in the playerLookDir vector, but only use the Y direction component
		testPlayerX -= playerLookDirY * playerStrafe;
		//testPlayerY += playerLookDirX * playerStrafe;

		testCellX = Math.floor(testPlayerX / cellDisplayScale);
		testCellY = Math.floor(testPlayerY / cellDisplayScale);
		testCellIndex = testCellY * numOfCellsX + testCellX;

		if (testCellX > -1 && testCellX < numOfCellsX && testCellY > -1 && testCellY < numOfCellsY &&
			cellsArray[testCellIndex] == 0)
		{ // success with moving in X direction only
			playerX = testPlayerX;
			//playerY = testPlayerY;
		}

		// reset testPlayer X and Y
		testPlayerX = playerX;
		testPlayerY = playerY;
		// trick to moving perpendicular to playerLookDir vector (sideways strafe)
		// move player's Y position in the playerLookDir vector, but only use the X direction component
		//testPlayerX -= playerLookDirY * playerStrafe;
		testPlayerY += playerLookDirX * playerStrafe;

		testCellX = Math.floor(testPlayerX / cellDisplayScale);
		testCellY = Math.floor(testPlayerY / cellDisplayScale);
		testCellIndex = testCellY * numOfCellsX + testCellX;

		if (testCellX > -1 && testCellX < numOfCellsX && testCellY > -1 && testCellY < numOfCellsY &&
			cellsArray[testCellIndex] == 0)
		{ // success with moving in Y direction only
			//playerX = testPlayerX;
			playerY = testPlayerY;
		}

	} // end if (KeyA_isPressed || KeyD_isPressed)

	// line of sight
	context.setLineDash([5, 10]);
	context.lineWidth = 1;
	context.strokeStyle = "cyan";
	context.moveTo(playerX, playerY);
	context.lineTo(mouseX, mouseY);
	context.stroke();
	// return line style to solid
	context.setLineDash([]);

	
	// player position
	context.lineWidth = 3;
	context.strokeStyle = "yellow";
	context.beginPath();
	context.arc(playerX, playerY, 10, playerLookAngle + 0.5, playerLookAngle - 0.5);
	context.stroke();
	
	// mouse lookAt target
	context.lineWidth = 2;
	context.strokeStyle = "blue";
	context.beginPath();
	context.arc(mouseX, mouseY, 10, 0, 2 * Math.PI);
	context.stroke();

	// ray intersection point (if any)
	if (tileFound)
	{
		context.lineWidth = 2;
		context.strokeStyle = "magenta";
		context.beginPath();
		context.arc(intersectionX, intersectionY, 10, 0, 2 * Math.PI);
		context.stroke();
	}

	

	debugInfo.innerHTML =
		"CellX: " + cellX + " CellY: " + cellY + " MouseCellIndex: " + mouseCellIndex + "<br>" +
		"PlayerCellX: " + playerCellX + " PlayerCellY: " + playerCellY + " PlayerCellIndex: " + playerCellIndex + "<br>" +
		"PlayerLookDirX: " + playerLookDirX.toPrecision(1) + " PlayerLookDirY: " + playerLookDirY.toPrecision(1) +
		" PlayerLookAngle: " + playerLookAngle.toPrecision(2);

	requestAnimationFrame(draw);
} // end function draw()

draw();
