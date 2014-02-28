
Modeling and Solving Uncertain Problems
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In the above example we had a simple box uncertainty set defined at the same time as the creation of the uncertainties. Modeling more complex uncertainty sets however requires additional machinery. First of all, let us consider polyhedral uncertainty sets.

By default JuMP treats all uncertainties as belonging to a polyhedral uncertainty set. You can build these polyhedral uncertainty sets using ``addConstraint`` with expressions only containing numbers and uncertainties, e.g. ``addConstraint(m, (2.0*u1-2.0) + (4.0*u2-2.0) <= +1)``.

While this may be sufficient for many situations, a rich variety of uncertainty sets exist, as well as at least three approaches for solving problems that can be combined in a variety of ways: reformulation, cutting planes, and sampling from the uncertainty set. In fact, we may want to use reformulate for one set of constraints with a polyhedral uncertainty set, and using cutting planes with an arbitrary convex uncertainty set for another set of constraints. To handle this complexity we introduce the concept of a **wrangler**.

A **wrangler** can be thought of as a "robustifying operator" that is responsible for taking one or many constraints stated in a raw, natural form and handling the task of ensuring that the solution is robust feasible with respect to the constraints under its command. (In fact, an explicit statement of the constraint is not absolutely required if it is not representable with JuMP.) To make this more tangible consider one of the wranglers that is provided with JuMP, ``PolyhedralWrangler``. This wrangler is the default wrangler and is applied to all constraints with uncertainties in them unless replaced with another. It implements the full interface that wranglers must support (see ``wrangler.jl`` for the full specification):

* **registerConstraint** - Notifies the wrangler that it is responsible for this constraint, and passes any preferences provided via the solve command. Returns a dictionary where the keys are the symbols ``:Cut``, ``:Reform``, and ``:Sample`` and the values are true or false, where true indicates the wrangler has selected these operations for this constraint. (not currently used)

* **setup** - Gives wrangler time to do any setup it needs to do. Called after all constraints have been registered. Examples of work that could be done here include transforming the uncertainty set and generating a cutting plane model. May be called multiple times - this should be handled gracefully by the wrangler.

* **generateCut** - Called in the main loop every iteration for every uncertain constraint, returns number of cuts added. Wrangler has full control over action taken with each constraint, and can chose to generate no cuts. In the case of ``PolyhedralWrangler`` the current solution to the master problem is used to set the objective of the cutting plane problem, the cutting plane problem is solved, and a new constraint is added if it would affect the current solution (i.e. robust feasiblity is violated.)

* **generateReform** - Called before the main loop, adds anything it wants to the model. In the case of ``PolyhedralWrangler``, it will reformulate all constraints under its control using duality.

In the above example a "master" problem is mentioned. This is the internal representation of the model used by JuMP in the solve process, and is in fact a normal JuMP model. It is constructed by copying all variables from the ``RobustModel`` as well as any constraints that are certain. The new variables and constraints from reformulation and cutting planes are added to this problem.

To manually specify a wrangler for a constraint, add it as an argument to addConstraint, e.g. ``addConstraint(m, u*x + w*y <= 3, BertsimasSimWrangler())``

Wranglers Provided
^^^^^^^^^^^^^^^^^^

JuMP comes with some wranglers that can be used directly and as inspiration for
your own custom uncertainty sets and solution technqiues.

* **PolyhedralWrangler** - supports polyhedral uncertainty sets defined using
  bounds on the uncertainties and constraints added involving only uncertains.
  Supports reformulation, cutting plane, a sampling method.
* **BertsimasSimWrangler** - implements the robust constraint variant detailed
  in the Bertsimas and Sim (2004) paper. Uses the bound information to determine
  the amount each uncertain can vary, but limits the number of uncertainties
  that can vary from their nominal value using a parameter Gamma.