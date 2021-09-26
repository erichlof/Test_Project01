let canvas, ctx;

function drawShape()
{

        // get the canvas element using the DOM
        canvas = document.getElementById('mycanvas');

        // use getContext to use the canvas for drawing
        ctx = canvas.getContext('2d');

        for (let i = 0; i < 10; i++)
        {
                ctx.lineWidth = 1 + i;
                ctx.beginPath();
                ctx.moveTo(5 + i * 14, 5);
                ctx.lineTo(5 + i * 14, 140);
                ctx.stroke();
        }

}

drawShape();