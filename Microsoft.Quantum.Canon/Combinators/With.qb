namespace Microsoft.Quantum.Canon {
	open Microsoft.Quantum.Primitive;

    operation With(outerOperation : (Qubit[] => ():Adjoint), innerOperation : (Qubit[] => ()), target : Qubit[])  : ()
    {
        body {  
            outerOperation(target);
            innerOperation(target);
            (Adjoint(outerOperation))(target);
        }
    }



    operation WithA(outerOperation : (Qubit[] => ():Adjoint), innerOperation : (Qubit[] => ():Adjoint), target : Qubit[])  : ()
    {
        body {  
            outerOperation(target);
            innerOperation(target);
            (Adjoint(outerOperation))(target);
        }

        adjoint auto
    }


    operation WithC(outerOperation : (Qubit[] => ():Adjoint), innerOperation : (Qubit[] => ():Controlled), target : Qubit[])  : ()
    {
        body {  
            outerOperation(target);
            innerOperation(target);
            (Adjoint(outerOperation))(target);
        }

        controlled(controlRegister) {
            outerOperation(target);
            (Controlled(innerOperation))(controlRegister, target);
            (Adjoint(outerOperation))(target);
        }
    }

    operation WithCA(outerOperation : (Qubit[] => ():Adjoint), innerOperation : (Qubit[] => ():Adjoint,Controlled), target : Qubit[])  : ()
    {
        body {  
            outerOperation(target);
            innerOperation(target);
            (Adjoint(outerOperation))(target);
        }

        adjoint auto
        controlled(controlRegister) {
            outerOperation(target);
            (Controlled(innerOperation))(controlRegister, target);
            (Adjoint(outerOperation))(target);
        }
        controlled adjoint auto
    }

    operation With1(outerOperation : (Qubit => ():Adjoint), innerOperation : (Qubit => ()), target : Qubit)  : ()
    {
        body {  
            outerOperation(target);
            innerOperation(target);
            (Adjoint(outerOperation))(target);
        }
    }

    operation With1C(outerOperation : (Qubit => ():Adjoint), innerOperation : (Qubit => ():Controlled), target : Qubit)  : ()
    {
        body {  
            With1(outerOperation, innerOperation, target);
        }

        controlled(controlRegister) {
            With1(outerOperation, (Controlled innerOperation)(controlRegister, _), target);
        }
    }

}