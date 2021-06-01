#if 1
#ifdef HSIMT
hsimt: block
             real(real64) :: kappa,x
                kappa = 1 - cfl
                x = 0.25_real64*( kappa - 1._real64/3._real64/kappa )
                limiter = (0.5_real64+x) + (0.5_real64-x)*ratio ! can become negative!!!
                limiter = MAX( 0.0_real64 , limiter )
                limiter = MIN( 2._real64*ratio , limiter , 2.0_real64 )
end block hsimt
#endif

#ifdef MUSCL
muscl: block
                limiter = MIN( 2*ratio , 0.5*(1+ratio) , 2.0_real64 )
end block muscl
#endif

#ifdef P2_PDM
ps_pdm: block
             real(real64) :: x
                x = 1._real64/6._real64*(1._real64-2._real64*cfl)
                limiter = (0.5_real64+x) + (0.5_real64-x)*ratio
                limiter = MIN( 2._real64*ratio/(cfl+1.e-10_real64) , limiter , 2._real64/(1._real64-cfl) )
end block ps_pdm
#endif

#ifdef SPLMAX13
splmax13: block
                limiter = MIN(2._real64*ratio,1._real64/3._real64*MAX( 1._real64+2._real64*ratio,2._real64+ratio ),2._real64)
end block splmax13
#endif

#ifdef SUPERBEE
superbee: block
                limiter = MAX( MIN( 2._real64*ratio , 1.0_real64 ) , MIN( ratio , 2.0_real64 ) )
end block superbee
#endif

#ifdef UPSTREAM
upstream: block
                limiter = 0._real64
end block upstream
#endif

#else
select case(scheme)
            case (1)
               hsimt: block
               real(real64) :: kappa,x
               kappa = 1 - cfl
               x = 0.25_real64*( kappa - 1._real64/3._real64/kappa )
               limiter = (0.5_real64+x) + (0.5_real64-x)*ratio ! can become negative!!!
               limiter = MAX( 0.0_real64 , limiter )
               limiter = MIN( 2._real64*ratio , limiter , 2.0_real64 )
               end block hsimt

            case (2)
               muscl: block
                limiter = MIN( 2*ratio , 0.5*(1+ratio) , 2.0_real64 )
               end block muscl

            case (3)
               p2_pdm: block
                real(real64) :: x
                x = 1._real64/6._real64*(1._real64-2._real64*cfl)
                limiter = (0.5_real64+x) + (0.5_real64-x)*ratio
                limiter = MIN( 2._real64*ratio/(cfl+1.e-10_real64) , limiter , 2._real64/(1._real64-cfl) )
               end block p2_pdm

            case (4)
               splmax13: block
                limiter = MIN( 2._real64*ratio , 1._real64/3._real64*MAX( 1._real64+2._real64*ratio , 2._real64+ratio ) , 2._real64 )
               end block splmax13

            case (5)
               superbee: block
                limiter = MAX( MIN( 2._real64*ratio , 1.0_real64 ) , MIN( ratio , 2.0_real64 ) )
               end block superbee

            case (6)
               upstream: block
                limiter = 0._real64
               end block upstream

            end select
#endif
